#include "raylib.h"
#include <cmath>
#include <random>

#if defined(PLATFORM_WEB)
    #include <emscripten/emscripten.h>
#endif

enum I_Method {
    EULER,
    RUNGEKUTTA,
    LEAPFROG
};
enum P_Distribution {
    ZERO,
    UNIFORM,
    BINOMIAL
};
enum Plot {
    NO_PLOT,
    KINETIC,
    POTENTIAL,
    ALL
};
enum Screen {
    START,
    SIMULATION
};

struct histogram {
    float rec;
    histogram *old;
};
struct block {
    int index = -1;
    long double x, v, v_prev, a;
    float e, k_p, k_c, friction;
    block *prev, *next;
};

int SCREEN_WIDTH  = 800;
int SCREEN_HEIGHT = 450;

// Parameters
const float k_p = 1.; // Molles amb el sostre
const float k_c = .4; // Molles entre blocs
const float friction_d = 10.;
const float v_e = 1.;
const float v_epsilon = 1e-3;
const float dt = 1e-5;
const int N_BLOCKS = 10;
const I_Method integrate = LEAPFROG;
const P_Distribution prob = BINOMIAL;
const Plot plot = ALL;
//const bool wrapping = false;

float max_rec = 100.;
histogram* kineticHist;
histogram* potentialHist;

block  *startBlock, *selected, *endBlock;

// Drawing settings
const int block_width = 3; 
float meter = 10.; // 1m = 10px
float x_displace = 0.;
float max_x;
bool follow = false;
bool paused = false;
Screen screen = START;

// Random engines
const int random_steps = 20;
std::default_random_engine generator;
std::uniform_int_distribution<int> uniform(0, random_steps);
std::binomial_distribution<int> binomial(random_steps, .5);

histogram* clear_histogram(histogram *start) {
    if (start != nullptr) {
        histogram *current = start->old;
        while (current != nullptr) {
            histogram *tmp = current;
            current = current->old;
            delete tmp;
        }
        delete start;
    }

    start = new histogram;
    start->rec = 0.;
    start->old = nullptr;
    return start;
}

void compute_starting_blocks() {
    if (plot == KINETIC or plot == ALL) kineticHist = clear_histogram(kineticHist);
    if (plot == POTENTIAL or plot == ALL) potentialHist = clear_histogram(potentialHist);

    if (startBlock != nullptr) {
        block* a = startBlock->next;
        while (a != nullptr) {
            block* tmp = a;
            a = a->next;
            delete tmp;
        }
        delete startBlock;
    }
    
    startBlock = new block;
    startBlock->prev = startBlock->next = nullptr;

    startBlock->index = 1;
    startBlock->x = startBlock->v = startBlock->a = 0.;
    startBlock->e = 2*block_width;
    startBlock->k_c = k_c;
    startBlock->k_p = k_p;
    startBlock->friction = friction_d;
    if (prob == UNIFORM) 
        startBlock->friction *= 1 + uniform(generator)/(float)random_steps;
    else if (prob == BINOMIAL) 
        startBlock->friction *= 1 + binomial(generator)/(float)random_steps;

    block *a = startBlock;
    for (int i = 1; i < N_BLOCKS; ++i) {
        block* b = new block;
        b->index = i + 1;
        b->x = 3*block_width*i;
        b->v = b->a = 0.;
        b->e = block_width*(3*i + 2);
        b->k_c = k_c;
        b->k_p = k_p;
        b->friction = friction_d;
        if (prob == UNIFORM) 
            b->friction *= 1 + uniform(generator)/(float)random_steps;
        else if (prob == BINOMIAL) 
            b->friction *= 1 + binomial(generator)/(float)random_steps;
        b->prev = a;
        b->next = nullptr;
        a = a->next = b;
    }
    max_x = a->x;

    selected = nullptr;
}

float get_gradient(float x, float min, float cutoff) {
    if (x < min) return 1 - (1-cutoff)/min*x;
    return cutoff;
}

float max(float a, float b) {
    if (a >= b) return a;
    else return b;
}

float mod(float a) {
    if (a < 0) return -a;
    return a;
}

int sign(float a) {
    if (a < 0) return -1;
    return 1;
}

Vector2 get_relative_mouse_pos() {
    Vector2 tmp = GetMousePosition();
    tmp.x += x_displace*meter;
    return tmp;
}

float compute_accel(const block& b) {
    float a = 0;
    if (b.next != nullptr) a = b.k_c*(b.next->x - b.x - block_width);
    if (b.prev != nullptr) a = b.prev->k_c*(b.prev->x - b.x - block_width);
    auto k_e_x = b.k_p*(b.e - b.x);
    if (mod(b.v) <= v_epsilon) {
        a = max(a + k_e_x - sign(b.v)*b.friction, 0);
    }
    else a += k_e_x - sign(b.v)*friction_d;
    return a;
}

float compute_vel(float dt, const block& b) {
    block tmp(b);
    float k_1 = b.a;
    float half_dt = dt/2;
    float v_1 = b.v + k_1*half_dt;
    tmp.x = b.x + v_1*half_dt;
    tmp.v = v_1;
    tmp.e = b.e + v_e*half_dt;
    float k_2 = compute_accel(tmp);
    float v_2 = b.v + k_2*half_dt;
    tmp.x = b.x + v_2*half_dt;
    tmp.v = v_2;
    float k_3 = compute_accel(tmp);
    float v_3 = b.v + k_3*dt;
    tmp.e = b.e + v_e*dt;
    float k_4 = compute_accel(tmp);
    return b.v + (k_1 + 2*k_2 + 2*k_3 + k_4)*dt/6;
}

void update_block(float dt, block &b) {
    if (integrate == LEAPFROG) {
        b.a = b.k_p*(b.e - b.x) - sign(b.v)*friction_d;
        if (b.next != nullptr) b.a += b.k_c*(b.next->x - b.x - block_width);
        if (b.prev != nullptr) b.a -= b.prev->k_c*(b.x - b.prev->x - block_width);

        b.v_prev = b.v;
        b.v = b.v_prev + dt*b.a;
        if (b.v*b.v_prev < 0 or mod(b.v) < v_epsilon) {
            if (mod(b.a) < b.friction) b.a = b.v = 0;
        }

        b.x += dt*b.v;
        b.e += v_e*dt;
    }
    else if (integrate == RUNGEKUTTA) {
        b.a = compute_accel(b);
        b.v = compute_vel(dt, b);
        b.x += b.v*dt;
        b.e += v_e*dt;
    } 
    else if (integrate == EULER){
        auto k_e_x = k_p*(b.e-b.x);
        if (b.v <= v_epsilon) {
            b.a += max(k_e_x - sign(b.v)*b.friction, 0);
        }
        else {
            b.a += k_e_x - sign(b.v)*friction_d;
        }
        b.v += b.a*dt;
        b.x += b.v*dt;
        b.e += v_e*dt;
    }
}

float spring_gradient(float x) {
    if (x < 0) x = -x;
    return (1/(x+1/0.3))*30+0.3;
}

void draw_block(const block &b, float b_width) {
    if (b.next != nullptr) {
        Vector2 a_pos = { (float)(b.x - x_displace)*meter + b_width, SCREEN_HEIGHT/2.f};
        Vector2 b_pos = { (float)(b.next->x - x_displace)*meter , SCREEN_HEIGHT/2.f};
        DrawLineEx(a_pos, b_pos, spring_gradient((float)b.next->x - (float)b.x - block_width)*(k_c/k_p)*0.3*meter, GRAY);
    }

    Vector2 square_pos = {(float)(b.x-x_displace)*meter, float((SCREEN_HEIGHT-b_width)/2)};

    Vector2 left_stick_dim = {meter/2, 6*meter};
    Vector2 left_stick_pos = {square_pos.x + b_width/2 - left_stick_dim.x/2, square_pos.y - left_stick_dim.y};

    Vector2 right_stick_dim(left_stick_dim);
    Vector2 right_stick_pos = {(float)(b.e-x_displace)*meter - right_stick_dim.x/2, square_pos.y - 10*meter};

    Color c;
    if (mod(b.v) > v_epsilon) c = GRAY;
    else c = RED;

    DrawRectangle(square_pos.x, square_pos.y, b_width, b_width, c);
    if (&b == selected) {
        DrawRectangleLinesEx((Rectangle){square_pos.x, square_pos.y, b_width, b_width}, 0.2*meter, BLUE);
    }

    Vector2 spring_start = {left_stick_pos.x + meter/2, left_stick_pos.y + 1*meter};
    Vector2 spring_end = {right_stick_pos.x, spring_start.y};

    DrawLineEx(spring_start, spring_end, spring_gradient((double)b.x+block_width/2.-(double)b.e)*0.3*meter, GRAY);
    DrawRectangleV(left_stick_pos, left_stick_dim, LIGHTGRAY);
    DrawRectangleV(right_stick_pos, right_stick_dim, LIGHTGRAY);
    
    if (b.index != -1) {
        DrawText(TextFormat("%i", b.index), square_pos.x+0.1*meter, square_pos.y, 2*meter, RAYWHITE);
    }

}

void draw_graph(histogram *start, Color color) {
    float height = (SCREEN_HEIGHT - block_width*meter)/2;
    if (start->rec > max_rec) max_rec = start->rec;

    const char* t = TextFormat("%.1f", fmax(start->rec, 0));
    // Marker text
    DrawText(t, SCREEN_WIDTH - MeasureText(t, 20), SCREEN_HEIGHT - start->rec/max_rec*height - 20, 20, color);
    DrawLineEx((Vector2){SCREEN_WIDTH - 20.f, SCREEN_HEIGHT - start->rec/max_rec*height}, (Vector2){(float)SCREEN_WIDTH, SCREEN_HEIGHT - start->rec/max_rec*height}, 3, color);

    int i = SCREEN_WIDTH;
    histogram *a = start;
    histogram *b = start->old;
    while (b != nullptr and i > 0) {
        DrawLine(i-1, SCREEN_HEIGHT - b->rec/max_rec*height, i, SCREEN_HEIGHT - a->rec/max_rec*height, color);
        a = b;
        b = b->old;
        --i;
    }

    a->old = nullptr;
    while (b != nullptr) {
        histogram *tmp = b;
        b = b->old;
        delete tmp;
    }
};

void draw_ui() {
    DrawFPS(0, 0);

    const char* z = TextFormat("Zoom: x%.3f", meter/10);
    DrawText(z, SCREEN_WIDTH-MeasureText(z, 20), 0, 20, GRAY);
    if (follow) {
        const char* t = "Follow";
        DrawText(t, SCREEN_WIDTH - MeasureText(z, 20) - MeasureText(t, 20) - 10, 0, 20, RED);
    }
    if (paused) {
        Vector2 pauseSize = {20, 70};
        Vector2 pausePos = {(SCREEN_WIDTH - pauseSize.x)/2 - 10, (SCREEN_HEIGHT-pauseSize.y)/2};
        DrawRectangleV(pausePos, pauseSize, DARKGRAY);
        pausePos.x += pauseSize.x + 20;
        DrawRectangleV(pausePos, pauseSize, DARKGRAY);
    }

    if (selected != nullptr) {
        int u_s = MeasureText(TextFormat("Upper Spring: %.2f", selected->k_p), 20);
        int r_s = MeasureText(TextFormat("Right Spring: %.2f", selected->k_c), 20);
        int k_e = MeasureText(TextFormat("Kinetic Energy: %.2f", selected->v*selected->v/2), 20);
        int f_s = MeasureText(TextFormat("Static Friction: %.2f", selected->friction), 20);
        int f_d = MeasureText(TextFormat("Dynamic Friction: %.2f", friction_d), 20);
        DrawText(TextFormat("Block %i:", selected->index), 100, 0, 20, GRAY);
        DrawText(TextFormat("Upper Spring: %.2f", selected->k_p), 100, 30, 20, GRAY);
        DrawText(TextFormat("Right Spring: %.2f", selected->k_c), 100 + u_s + 20, 30, 20, GRAY);
        DrawText(TextFormat("Kinetic Energy: %.2Lf", selected->v*selected->v/2), 100 + u_s + r_s + 40, 30, 20, GRAY);
        DrawText(TextFormat("Static Friction: %.2f", selected->friction), 100, 60, 20, GRAY);
        DrawText(TextFormat("Dynamic Friction: %.2f", friction_d), 100 + f_s + 20, 60, 20, GRAY);
    }
}

void UpdateDrawFrame();

int main()
{
    InitWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Earthquake");

    compute_starting_blocks();

#if defined(PLATFORM_WEB)
    emscripten_set_main_loop(UpdateDrawFrame, 0, 1);
#else
    SetTargetFPS(60);
    while (!WindowShouldClose()) UpdateDrawFrame();
#endif

    CloseWindow();

    return 0;
}

void UpdateDrawFrame(void)
{
    SCREEN_WIDTH = GetScreenWidth();
    SCREEN_HEIGHT = GetScreenHeight();

    if (screen == START) {
        const char* sc = "Source Code https://github.com/bibanez/earthquake-sim";
        Vector2 start = {50, 50};
        if (GetKeyPressed() != 0) screen = SIMULATION;
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
            Vector2 pos = GetMousePosition();
            if (pos.x >= start.x and pos.x <= start.x + MeasureText(sc, 20)
                    and pos.y >= start.y + 300 and pos.y <= start.y + 320) {
                OpenURL("https://github.com/bibanez/earthquake-sim");
            }
        }

        BeginDrawing();

        ClearBackground(RAYWHITE);

        DrawText("The purpose of this simulation is to play", start.x, start.y, 20, GRAY);
        DrawText("with a simple dynamic model that shows some", start.x, start.y + 30, 20, GRAY);
        DrawText("of the characteristics of earthquakes.", start.x, start.y + 60, 20, GRAY);
        DrawText("Keybinds: [f] -> follow [Up/Down Arrows] -> Zoom", start.x, start.y + 100, 20, GRAY);
        DrawText("[Left/Right Arrows] -> Move horizontally", start.x, start.y + 130, 20, GRAY);
        DrawText("[b] -> back to menu  [r] -> reset simulation", start.x, start.y + 160, 20, GRAY);

        DrawText("Legend:", start.x, start.y + 270, 20, DARKGRAY);
        DrawText("Potential Energy", start.x + 100, start.y + 270, 20, DARKGREEN);
        DrawText("Kinetic Energy", start.x + 300, start.y + 270, 20, DARKBLUE);

        DrawText("Press any key to start", start.x, start.y + 240, 20, GRAY);
        DrawText(sc, start.x, start.y + 320, 20, BLUE);


        const char* name = "Made by bibanez (Bernat Ibáñez)";
        DrawText(name, SCREEN_WIDTH - MeasureText(name, 20) - 10, SCREEN_HEIGHT - 30, 20, GRAY);
        EndDrawing();
    }
    else { //  if (screen == SIMULATION)
        char c = GetCharPressed();
        if (c == 'r') {
            x_displace = 0;
            compute_starting_blocks();
            meter = 10.;
        }
        else if (c == '+' or IsKeyDown(KEY_UP)) meter += 0.1;
        else if ((c == '-' or IsKeyDown(KEY_DOWN)) and meter > 0.1) meter -= 0.1;
        else if (c == 'f') follow = !follow;
        else if (c == 'b') {
            screen = START;
            follow = paused = false;
            compute_starting_blocks();
        }
        if (c == 'p' or IsKeyPressed(KEY_SPACE)) paused = !paused;
        
        if (IsKeyDown(KEY_LEFT)) {
            follow = false;
            x_displace -= 1;
        }
        if (IsKeyDown(KEY_RIGHT)) {
            follow = false;
            x_displace += 1;
        }

        auto b_width = block_width*meter;

        Vector2 mouse = get_relative_mouse_pos();
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
            if (mouse.y <= (SCREEN_HEIGHT + b_width)/2 and mouse.y >= (SCREEN_HEIGHT - b_width)/2) {
                block* a = startBlock;
                while (a != nullptr and (mouse.x < a->x*meter or mouse.x > a->x*meter + b_width)) a = a->next;
                selected = a;
            }
            else selected = nullptr;
        }
        

        if (not paused) {
            int iter = GetFrameTime()/dt;
            for (int i = 0; i < iter; ++i) {
                block *b = startBlock;
                while (b != nullptr) {
                    if (b->x > max_x) max_x = b->x;
                    update_block(dt, *b);
                    b = b->next;
                }
            }

            if (plot == KINETIC or plot == ALL) {
                histogram *tmp = kineticHist;
                kineticHist = new histogram;
                kineticHist->rec = 0.;
                kineticHist->old = tmp;

                block *b = startBlock;
                while (b != nullptr) {
                    if (LEAPFROG) {
                        float v = (b->v + b->v_prev)/2;
                        kineticHist->rec += v/2;
                    }
                    else kineticHist->rec += b->v*b->v/2;

                    b = b->next;
                }
            } 
            if (plot == POTENTIAL or plot == ALL) {
                histogram *tmp = potentialHist;
                potentialHist = new histogram;
                potentialHist->rec = 0.;
                potentialHist->old = tmp;

                block *b = startBlock;
                while (b != nullptr) {
                    float x = (b->e - b->x);
                    potentialHist->rec += b->k_p*x*x/2;
                    b = b->next;
                }
            }
        }
        if (follow) x_displace = max_x + block_width - SCREEN_WIDTH/meter;
        

        BeginDrawing();

        ClearBackground(RAYWHITE);
       
        DrawRectangle(0, (SCREEN_HEIGHT + b_width)/2, SCREEN_WIDTH, SCREEN_HEIGHT, LIGHTGRAY);
        DrawRectangle(0, (SCREEN_HEIGHT - b_width)/2 - 10*meter, SCREEN_WIDTH, 1*meter, LIGHTGRAY);

        block *b = startBlock;
        while (b != nullptr) {
            draw_block(*b, b_width);
            b = b->next;
        }

        if (plot == KINETIC or plot == ALL) 
            draw_graph(kineticHist, DARKBLUE);
        if (plot == POTENTIAL or plot == ALL) 
            draw_graph(potentialHist, DARKGREEN);

        draw_ui();

        EndDrawing();
    }
    }
