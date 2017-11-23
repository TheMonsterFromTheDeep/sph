#include <GLFW/glfw3.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>

#define REF_DENSITY 100
#define BULK_MODULUS 2
#define VISCOSITY 0.1f
#define GRAV_X 0
#define GRAV_Y -9.8f
#define PARTICLE_SIZE 5e-2

float mass;

typedef struct {
    float x, y;
    float vx, vy;
    float ax, ay;
    float density;
} particle;

void get_densities(particle *list, size_t n) {
    float sz2 = PARTICLE_SIZE * PARTICLE_SIZE;
    float C = (4 * mass) / (M_PI * (sz2 * sz2) * (sz2 * sz2));
    
    for(size_t i = 0; i < n; ++i) {
        list[i].density = 0;
    }
    
    for(size_t i = 0; i < n; ++i) {
        /* Contribution from itself */
        list[i].density += (4 * mass) / (M_PI * sz2);
        for(size_t j = i + 1; j < n; ++j) {
            float dx = list[i].x - list[j].x;
            float dy = list[i].y - list[j].y;
            float dist2 = dx * dx + dy * dy;
            float dif = sz2 - dist2;
            if(dif > 0) {
                float contrib = C * dif * dif * dif;
                list[i].density += contrib;
                list[j].density += contrib;
            }
        }
    }
}

void get_accelerations(particle *list, size_t n) {
    float sz2 = PARTICLE_SIZE * PARTICLE_SIZE;
    float C = mass / (M_PI * sz2 * sz2);
    float k15 = 15 * BULK_MODULUS;
    float mu40 = -40 * VISCOSITY;
    
    for(size_t i = 0; i < n; ++i) {
        /* Begin acceleration calculation with gravitational force */
        list[i].ax = GRAV_X;
        list[i].ay = GRAV_Y;
    }
    
    for(size_t i = 0; i < n; ++i) {
        for(size_t j = i + 1; j < n; ++j) {
            float dx = list[i].x - list[j].x;
            float dy = list[i].y - list[j].y;
            float dist2 = dx * dx + dy * dy;

            if(dist2 < sz2) {
                float idens = list[i].density;
                float jdens = list[j].density;
                
                float norm = (sqrt(dist2) / PARTICLE_SIZE);
                float mnorm = 1 - norm;
                
                float cof = (C * mnorm) / (jdens * idens);
                float kterm = (cof * k15 * (idens + jdens - 2 * REF_DENSITY)
                                  * mnorm) / norm;
                float muterm = cof * mu40;
                
                float fx = kterm * dx + muterm * (list[i].vx - list[j].vx);
                float fy = kterm * dy + muterm * (list[i].vy - list[j].vy);
                
                list[i].ax += fx;
                list[i].ay += fy;
                list[j].ax -= fx;
                list[j].ay -= fy;
            }
        }
    }
}

#define DAMP 0.1f

void reflect_x(particle *p, float where) {
    //if(p->vx == 0) return;
    p->x = where;
    p->vx *= -DAMP;
}

void reflect_y(particle *p, float where) {
    //if(p->vy == 0) return;
    p->y = where;
    p->vy *= -DAMP;
}

#define XMIN 0.05f
#define XMAX 0.95f
#define YMIN 0.05f

void integrate_and_bound(particle *list, size_t n, float dt) {
    for(size_t i = 0; i < n; ++i) {
        list[i].vx += list[i].ax * dt;
        list[i].vy += list[i].ay * dt;
        list[i].x += list[i].vx * dt;
        list[i].y += list[i].vy * dt;
        
        if(list[i].x < XMIN) reflect_x(list + i, XMIN);
        if(list[i].x > XMAX) reflect_x(list + i, XMAX);
        if(list[i].y < YMIN) reflect_y(list + i, YMIN);
    }
    
    int good = 1;
    
    for(size_t i = 0; i < n; ++i) {
        if(list[i].x < XMIN) good = 0;
        if(list[i].x > XMAX) good = 0;
        if(list[i].y < YMIN) good = 0;
    }
}

void step(particle *list, size_t n, float dt) {
    get_densities(list, n);
    get_accelerations(list, n);
    integrate_and_bound(list, n, dt);
}

void draw(particle *p) {
    glBegin(GL_QUADS);
        glVertex2f(p->x - (PARTICLE_SIZE / 2), p->y);
        glVertex2f(p->x, p->y - (PARTICLE_SIZE / 2));
        glVertex2f(p->x + (PARTICLE_SIZE / 2), p->y);
        glVertex2f(p->x, p->y + (PARTICLE_SIZE / 2));
    glEnd();
}

int indicator(float x, float y) {
    x -= 0.5f;
    y -= 0.5f;
    return x * x + y * y <= 0.4f * 0.4f;
}

void calculate_mass(particle *list, size_t n) {
    /* Calculate initial densities, pretending mass is 1 */
    mass = 1;
    get_densities(list, n);
    
    float dens = 0;
    float dens2 = 0;
    for(size_t i = 0; i < n; ++i) {
        dens += list[i].density;
        dens2 += list[i].density * list[i].density;
    }
    
    mass = (REF_DENSITY * dens) / dens2;
}

void populate_list(particle **list_out, size_t *size) {
    size_t count = 0;
    
    float delta = PARTICLE_SIZE / 1.3f;
    
    for(float x = 0; x < 1; x += delta) {
        for(float y = 0; y < 1; y += delta) {
            count += indicator(x, y);
        }
    }
    
    particle *list = malloc(sizeof(particle) * count);
    
    for(size_t i = 0; i < count; ++i) {
        list[i].x = list[i].y = list[i].ax = list[i].ay = list[i].vx = list[i].vy = 0;
    }
    
    size_t i = 0;
    
    for(float x = 0; x < 1; x += delta) {
        for(float y = 0; y < 1; y += delta) {
            if(indicator(x, y)) {
                list[i].x = x;
                list[i].y = y;
                ++i;
            }
        }
    }
    
    calculate_mass(list, count);
    
    *size = count;
    *list_out = list;
}

#define GRID_SIZE 300

float density_grid[GRID_SIZE][GRID_SIZE];

void grid_particle(particle *p) {
    size_t mx = (size_t)(p->x * GRID_SIZE);
    size_t my = (size_t)(p->y * GRID_SIZE);
    int extent = (int)(PARTICLE_SIZE * GRID_SIZE / 2);
    float scaf = p->density / REF_DENSITY;
    
    for(int x = mx - extent; x < mx + extent; ++x) {
        if(x < 0) continue;
        if(x >= GRID_SIZE) break;
        for(int y = my - extent; y < my + extent; ++y) {
            if(y < 0) continue;
            if(y >= GRID_SIZE) break;
            int dx = x - mx;
            int dy = y - my;
            int r = dx * dx + dy * dy;
            if(r <= extent * extent) {
                float dens = 1 - (float)r / (extent * extent);
                density_grid[x][y] += dens * scaf;
            }
        }
    }
}

void fill_grid(particle *list, size_t n) {
    for(size_t x = 0; x < GRID_SIZE; ++x) {
        for(size_t y = 0; y < GRID_SIZE; ++y) {
            density_grid[x][y] = 0;
        }
    }
    
    for(size_t i = 0; i < n; ++i) {
        grid_particle(list + i);
    }
    
    /*for(size_t x = 0; x < 100; ++x) {
        for(size_t y = 0; y < 100; ++y) {
            density_grid[x][y] /= 5;
        }
    }*/
}

void checker_grid() {
    for(size_t x = 0; x < GRID_SIZE; ++x) {
        for(size_t y = 0; y < GRID_SIZE; ++y) {
            density_grid[x][y] = x % 2;
        }
    }
}

void render_grid() {
    for(size_t x = 0; x < GRID_SIZE; ++x) {
        for(size_t y = 0; y < GRID_SIZE; ++y) {
            
            float cx = (float)x / GRID_SIZE;
            float cy = (float)y / GRID_SIZE;
            float p = (float)1 / GRID_SIZE;
            
            float c = density_grid[x][y];
            if(c < 0) c = 0;
            if(c > 1) c = 1;
            
            glColor3f(c, c, c);
            
            glBegin(GL_QUADS);
                glVertex2f(cx, cy);
                glVertex2f(cx, cy + p);
                glVertex2f(cx + p, cy + p);
                glVertex2f(cx + p, cy);
            glEnd();
        }
    }
}

int main(void) {
    GLFWwindow* window;

    if (!glfwInit())
        return -1;

    window = glfwCreateWindow(640, 640, "Fluid Simulation Test", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glOrtho(0, 1, 0, 1, -20, 20);
    
    size_t particle_count;
    particle *fluid;
    
    populate_list(&fluid, &particle_count);
    
    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);
        
        render_grid();
        
        glfwSwapBuffers(window);
        glfwPollEvents();
        
        step(fluid, particle_count, 0.01f);
        fill_grid(fluid, particle_count);
    }
    
    free(fluid);

    glfwTerminate();
    return 0;
}