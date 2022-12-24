#define DEBUG

#include <cstdlib>
#include <string>
#include <random>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <thread>
#include "Resources\EasyX\include\graphics.h"
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#pragma comment(lib, "winmm.lib")

#ifdef DEBUG
#include <iostream>
#endif

/** Global variable declarations and initializations */
WCHAR working_path[MAX_PATH];
DWORD* graph_buffer /* display buffer of the graph */;
long long page;
short current_bgm;
#define GRAPH_SIZE 1000 /* (in pixels) avoid changing it */
#define SHUFFLE_TIMES 70000 /* for TRIANGLE_MODE, use a number close to the square of the number of blocks */
/** -----------------------------------------------  */

/** Mode and page definitions */
/*                        x      */
#define START_PAGE      0b100000
#define EXIT_PAGE       0b100001
/*                         x     */
#define SQR_TRI_FILTER  0b010000 /* bitfilter */
#define SQUARE_MODE     0b000000
#define TRIANGLE_MODE   0b010000
/*                          x    */
#define HUE_IMG_FILTER  0b001000 /* bitfilter */
#define HUE_MODE        0b000000
#define IMAGE_MODE      0b001000
/*                           x   */
#define RESUMED_FILTER  0b000100 /* bitfilter */
#define RESUMED         0b000100
/*                            c_ (page type) */
#define PAGE_TYPE       0b000011 /* bitfilter */
#define GAME_PAGE       0b000000
#define PAUSED_PAGE     0b000001
#define RESULT_PAGE     0b000010
#define RANK_PAGE       0b000011

//TODO #define SETTINGS_PAGE   0b000001
/** ------------------------ */

class Triangle {
public:
    std::pair<long long, long long> vertices[3] /* (Ax, Ay), (Bx, By), (Cx, Cy) */;
    short type /* 0 for upward and 1 for downward (take the primitive triangle as upward) */, boundary /* a number between 0 and 7 (inclusive), treated as a binary number */;
    long long target /* a number between 0 and custom::raise_to_power(4, ?) - 1 (inclusive), to be changed when swapping */, direction /* a number between 0 and 2 (inclusive), to be changed when swapping */;
};

class Square {
public:
    long long begin_x, begin_y, size;
    long long target /* a number between 0 and N * N - 1 (inclusive), to be changed when swapping */;
};

class Hue {
public:
    int top_left_RGB[3], total_change_i[3] /* the total change in the i-direction */, total_change_j[3] /* the total change in the j-direction */;
};


namespace custom {
    long long rand(long long M) /* return an integer between 0 and M - 1 (inclusive) */ {
        std::random_device rd1 /* a seed for the random number engine */;
        std::mt19937 rd2(rd1()) /* standard mersenne_twister_engine seeded with rd1() */;
        // Use either rd1 or rd2 as the seed to feed the distribution below
        std::uniform_int_distribution<long long> distrib(0, M - 1);
        return distrib(rd2);
    }
    long long raise_to_power(long long a, long long b /* non-negative integer */) {
        if (b < 0)
            throw std::runtime_error("Negative exponent is not supported.");
        long long ans = 1;
        while (b)
            ans *= a, --b;
        return ans;
    }
    bool comp_x(const std::pair<long long, long long> A, const std::pair<long long, long long> B) {
        return (A.first < B.first) || (A.first == B.first && A.second < B.second);
    }
    bool comp_y(const std::pair<long long, long long> A, const std::pair<long long, long long> B) {
        return (A.second < B.second) || (A.second == B.second && A.first < B.first);
    }
    std::wstring full_path(const wchar_t* relative_path) {
        /** Make sure initialization() is first runned */
        std::wstring wstrpath(working_path);
        wstrpath += relative_path;
        return wstrpath;
    }
    std::string full_path(const char* relative_path) {
        /** Make sure initialization() is first runned */
        std::wstring wstrpath(working_path);
        std::string strpath(wstrpath.begin(), wstrpath.end());
        strpath += relative_path;
        return strpath;
    }
    bool is_in_rect(ExMessage* msg, int x, int y, int delta_x, int delta_y) {
        return x <= msg->x && msg->x < x + delta_x && y <= msg->y && msg->y < y + delta_y;
    }
    bool is_in_circ(ExMessage* msg, int x, int y, int radius) {
        return (msg->x - x) * (msg->x - x) + (msg->y - y) * (msg->y - y) < radius * radius;
    }
    long long consolas_width(long long height) {
        long double width = height * 0.493L /* width to height ratio of the "Consolas" font */;
        return (long long)(width + 0.5);
    }
}

namespace square {
    void swap(Square* A, Square* B) {
        long long tmp = A->target;
        A->target = B->target;
        B->target = tmp;
    }
}

namespace triangle {
    void swap(Triangle* A, Triangle* B, long long adjacent_type) {
        Triangle* tmp = new Triangle;
        *tmp = *A;
        A->target = B->target;
        A->direction = (6 - adjacent_type - B->direction) % 3;
        B->target = tmp->target;
        B->direction = (6 - adjacent_type - tmp->direction) % 3;
    }
    bool swap_if_adjacent(Triangle* A, Triangle* B) /* return true if the triangles are adjacent and swapping is successfully done */ {
        for (short k = 0; k != 3; ++k)
            if (A->vertices[(k + 1) % 3] == B->vertices[(k + 2) % 3] && A->vertices[(k + 2) % 3] == B->vertices[(k + 1) % 3]) {
                swap(A, B, k);
                return true;
            }
        return false;
    }
}

inline const int color_hue_BGR(const int i, const int j, const Hue* hue) {
    return BGR(RGB(hue->top_left_RGB[0] + i * hue->total_change_i[0] / GRAPH_SIZE + j * hue->total_change_j[0] / GRAPH_SIZE, hue->top_left_RGB[1] + i * hue->total_change_i[1] / GRAPH_SIZE + j * hue->total_change_j[1] / GRAPH_SIZE, hue->top_left_RGB[2] + i * hue->total_change_i[2] / GRAPH_SIZE + j * hue->total_change_j[2] / GRAPH_SIZE));
}

class Game {
private:
    short mode;
    long long N, empty_block, original_empty_block, game_move;
    std::chrono::milliseconds game_time;
    std::chrono::seconds game_time_limit /* to be assigned with a positive value */;
public:
    Square prim_sqr /* primitive (original non-splited) square */;
    Triangle prim_tri /* primitive (original non-splited) triangle */;
    Square* square_list /* array to be dynamically created */;
    Triangle** triangle_list /* array of arrays to be dynamically created */;
    long long* pixel_status /* useful in TRIANGLE_MODE */;
    Game(short mode, long long N /* positive integer */) {
        /** Initialize the game enviornment except that the size of prim_sqr or prim_tri is not given */
        if (N <= 0)
            throw std::runtime_error("N should be positive.");
        Game::N = N;
        Game::game_time_limit = std::chrono::seconds(0);
        Game::mode = mode;
        if (mode == TRIANGLE_MODE) {
            Game::triangle_list = new Triangle * [N];
            for (long long i = 0; i != N; ++i)
                Game::triangle_list[i] = new Triangle[custom::raise_to_power(4, i)];
            **Game::triangle_list = Game::prim_tri;
            pixel_status = new long long[GRAPH_SIZE * GRAPH_SIZE];
        }
        else if (Game::mode == SQUARE_MODE) {
            Game::square_list = new Square[N * N];
        }
        else {
            throw std::runtime_error("Not supported.");
        }
        /** ----------------------------------------------------------------------------------------- */
    }
    short get_mode() {
        return Game::mode;
    }
    long long get_N() {
        return Game::N;
    }
    long long get_empty(bool original = false) {
        return original? Game::original_empty_block: Game::empty_block;
    }
    std::chrono::milliseconds get_timer() {
        return Game::game_time;
    }
    void reset_timer() {
        Game::game_time = std::chrono::milliseconds(0);
    }
    void increment_timer(std::chrono::milliseconds t) {
        Game::game_time += t;
    }
    void finalize_timer() {
        if (Game::game_time_limit > std::chrono::milliseconds(0))
            Game::game_time = Game::game_time_limit;
        else
            throw std::runtime_error("Time limit not set.");
    }
    std::chrono::seconds get_time_limit() {
        return Game::game_time_limit;
    }
    void set_time_limit(std::chrono::seconds t) {
        Game::game_time_limit = t;
    }
    long long get_cnt() {
        return Game::game_move;
    }
    void reset_cnt() {
        Game::game_move = 0;
    }
    void increment_cnt() {
        Game::game_move++;
    }
    void set_prim_sqr(long long begin_x, long long begin_y, long long size) {
        Game::prim_sqr.begin_x = begin_x, Game::prim_sqr.begin_y = begin_y;
        Game::prim_sqr.size = size / Game::N * Game::N; //TODO maybe changed quite strangely when such permanant changes was applied for many times?
        if (!Game::prim_sqr.size)
            throw std::runtime_error("N is too large.");
    }
    void set_prim_tri(long long Ax, long long Ay, long long Bx, long long By, long long Cx, long long Cy) {
        Game::prim_tri.vertices[0] = std::pair<long long, long long>(Ax, Ay), Game::prim_tri.vertices[1] = std::pair<long long, long long>(Bx, By), Game::prim_tri.vertices[2] = std::pair<long long, long long>(Cx, Cy);
        Game::prim_tri.direction = 0, Game::prim_tri.type = 0, Game::prim_tri.boundary = 0b111; /* not to be changed */
    }
    void set_tri_recursion(long long layer, long long k, bool initialize) {
        if (initialize)
            Game::triangle_list[layer][k].target = k;
        if (layer == Game::N)
            return;
        std::pair<long long, long long> A = Game::triangle_list[layer][k].vertices[0], B = Game::triangle_list[layer][k].vertices[1], C = Game::triangle_list[layer][k].vertices[2];
        std::pair<long long, long long> F((A.first + B.first) >> 1, (A.second + B.second) >> 1), E((C.first + A.first) >> 1, (C.second + A.second) >> 1), D((B.first + C.first) >> 1, (B.second + C.second) >> 1);
        Game::triangle_list[layer + 1][4 * k].vertices[0] = A;
        Game::triangle_list[layer + 1][4 * k].vertices[1] = F;
        Game::triangle_list[layer + 1][4 * k].vertices[2] = E;
        if (initialize) {
            Game::triangle_list[layer + 1][4 * k].direction = Game::triangle_list[layer][k].direction;
            Game::triangle_list[layer + 1][4 * k].type = Game::triangle_list[layer][k].type;
            Game::triangle_list[layer + 1][4 * k].boundary = Game::triangle_list[layer][k].boundary & 0b110;
        }
        Game::set_tri_recursion(layer + 1, 4 * k, initialize);
        Game::triangle_list[layer + 1][4 * k + 1].vertices[0] = F;
        Game::triangle_list[layer + 1][4 * k + 1].vertices[1] = B;
        Game::triangle_list[layer + 1][4 * k + 1].vertices[2] = D;
        if (initialize) {
            Game::triangle_list[layer + 1][4 * k + 1].direction = Game::triangle_list[layer][k].direction;
            Game::triangle_list[layer + 1][4 * k + 1].type = Game::triangle_list[layer][k].type;
            Game::triangle_list[layer + 1][4 * k + 1].boundary = Game::triangle_list[layer][k].boundary & 0b101;
        }
        Game::set_tri_recursion(layer + 1, 4 * k + 1, initialize);
        Game::triangle_list[layer + 1][4 * k + 2].vertices[0] = E;
        Game::triangle_list[layer + 1][4 * k + 2].vertices[1] = D;
        Game::triangle_list[layer + 1][4 * k + 2].vertices[2] = C;
        if (initialize) {
            Game::triangle_list[layer + 1][4 * k + 2].direction = Game::triangle_list[layer][k].direction;
            Game::triangle_list[layer + 1][4 * k + 2].type = Game::triangle_list[layer][k].type;
            Game::triangle_list[layer + 1][4 * k + 2].boundary = Game::triangle_list[layer][k].boundary & 0b011;
        }
        Game::set_tri_recursion(layer + 1, 4 * k + 2, initialize);
        Game::triangle_list[layer + 1][4 * k + 3].vertices[0] = D;
        Game::triangle_list[layer + 1][4 * k + 3].vertices[1] = E;
        Game::triangle_list[layer + 1][4 * k + 3].vertices[2] = F;
        if (initialize) {
            Game::triangle_list[layer + 1][4 * k + 3].direction = Game::triangle_list[layer][k].direction;
            Game::triangle_list[layer + 1][4 * k + 3].type = (Game::triangle_list[layer][k].type + 1) % 2;
            Game::triangle_list[layer + 1][4 * k + 3].boundary = 0;
        }
        Game::set_tri_recursion(layer + 1, 4 * k + 3, initialize);
    }
    void set_tri(bool initialize) {
        Game::set_tri_recursion(Game::N - 1, 0, initialize);
    }
    void set_sqr(bool initialize) {
        for (long long i = 0; i != Game::N; ++i) {
            for (long long j = 0; j != Game::N; ++j) {
                Game::square_list[i * Game::N + j].begin_x = Game::prim_sqr.begin_x + j * (Game::prim_sqr.size / Game::N);
                Game::square_list[i * Game::N + j].begin_y = Game::prim_sqr.begin_y + i * (Game::prim_sqr.size / Game::N);
                Game::square_list[i * Game::N + j].size = Game::prim_sqr.size / Game::N /* maybe not used */;
                if (initialize)
                    Game::square_list[i * Game::N + j].target = i * Game::N + j;
            }
        }
    }
    bool check_solved() {
        bool winned = true;
        if (Game::mode == TRIANGLE_MODE) {
            for (long long i = 0; winned && i != custom::raise_to_power(4, Game::N - 1); ++i)
                winned = winned && Game::triangle_list[Game::N - 1][i].target == i;
        }
        else if (Game::mode == SQUARE_MODE)
            for (long long i = 0; winned && i != Game::N * Game::N; ++i)
                winned = winned && Game::square_list[i].target == i;
        return winned;
    }
    void move(long long delta_i, long long delta_j) {
        //TODO if (Game::mode == TRIANGLE_MODE) {
        Game::empty_block += delta_i * Game::N + delta_j;
    }
    void set_original_empty_block(long long empty) {
        Game::original_empty_block = empty;
    }
    void shuffle(long long shuffle_times, bool random_assign_original_empty_block = true) /* call this function after initialization */ {
        if (!Game::check_solved())
            throw std::runtime_error("Not initialized.");
        if (Game::mode == TRIANGLE_MODE) {
            if (random_assign_original_empty_block)
                Game::original_empty_block = custom::rand(custom::raise_to_power(4, Game::N - 1));
            Game::empty_block = Game::original_empty_block;
            for (long long k = 0; k != shuffle_times; ++k) /* the well behavior of this routine depends somehow on the precision of floating point arithmetic which are assumed*/ { //TODO Need to consider the case where no shuffle is done
                long long rand_dir = custom::rand(3);
                std::pair<long double, long double> empty_centroid((Game::triangle_list[Game::N - 1][Game::empty_block].vertices[0].first + Game::triangle_list[Game::N - 1][Game::empty_block].vertices[1].first + Game::triangle_list[Game::N - 1][Game::empty_block].vertices[2].first) / 3.0L, (Game::triangle_list[Game::N - 1][Game::empty_block].vertices[0].second + Game::triangle_list[Game::N - 1][Game::empty_block].vertices[1].second + Game::triangle_list[Game::N - 1][Game::empty_block].vertices[2].second) / 3.0F);
                std::pair<long long, long long> target_centroid(2.0L * empty_centroid.first - Game::triangle_list[Game::N - 1][empty_block].vertices[rand_dir].first - 0.5L, 2.0L * empty_centroid.second - Game::triangle_list[Game::N - 1][Game::empty_block].vertices[rand_dir].second - 0.5L);
                long long target;
                if (target_centroid.first >= 0 && target_centroid.first < GRAPH_SIZE && target_centroid.second >= 0 && target_centroid.second < GRAPH_SIZE && (target = Game::pixel_status[target_centroid.second * GRAPH_SIZE + target_centroid.first]) != -1) {
                    if (target == Game::empty_block /* I believe this should not happen */)
                        throw std::runtime_error("Should not equal. empty_block = " + std::to_string(Game::empty_block) + "; target = " + std::to_string(target) + ".\n");
                    else if (!triangle::swap_if_adjacent(*(Game::triangle_list + Game::N - 1) + Game::empty_block, *(Game::triangle_list + Game::N - 1) + target)) /* I believe the line below should not execute */
                        throw std::runtime_error("Not adjacent. empty_block = " + std::to_string(Game::empty_block) + "; target = " + std::to_string(target));
                    Game::empty_block = target;
                }
            }
        }
        else if (Game::mode == SQUARE_MODE) {
            if (random_assign_original_empty_block)
                Game::original_empty_block = custom::rand(Game::N * Game::N);
            Game::empty_block = Game::original_empty_block;
            for (long long k = 0; k != shuffle_times; ++k) { //TODO Need to consider the case where no shuffle is done
                long long rand_dir = custom::rand(4), delta_i, delta_j;
                switch (rand_dir) {
                case 0:
                    delta_i = 0, delta_j = 1;
                    break;
                case 1:
                    delta_i = 0, delta_j = -1;
                    break;
                case 2:
                    delta_i = 1, delta_j = 0;
                    break;
                case 3:
                    delta_i = -1, delta_j = 0;
                }
                if (delta_i + Game::empty_block / Game::N >= 0 && delta_i + Game::empty_block / Game::N != Game::N && delta_j + Game::empty_block % Game::N >= 0 && delta_j + Game::empty_block % Game::N != Game::N) {
                    square::swap(Game::square_list + Game::empty_block, Game::square_list + Game::empty_block + delta_i * Game::N + delta_j);
                    move(delta_i, delta_j);
                }
            }
        }
    }
};

namespace square {
    void fill_block_target(Game* game, const Hue* hue, long long block_i, long long block_j, long long target_block_i, long long target_block_j, bool empty) {
        for (long long i = 0; i != game->prim_sqr.size / game->get_N(); ++i)
            for (long long j = 0; j != game->prim_sqr.size / game->get_N(); ++j)
                graph_buffer[(game->prim_sqr.begin_y + block_i * (game->prim_sqr.size / game->get_N()) + i) * GRAPH_SIZE + (game->prim_sqr.begin_x + block_j * (game->prim_sqr.size / game->get_N()) + j)] = empty ? WHITE : color_hue_BGR(game->prim_sqr.begin_y + target_block_i * (game->prim_sqr.size / game->get_N()) + i, game->prim_sqr.begin_x + target_block_j * (game->prim_sqr.size / game->get_N()) + j, hue);
    }
    void fill_block(Game* game, const Hue* hue, long long block_i, long long block_j, bool empty) {
        fill_block_target(game, hue, block_i, block_j, (game->square_list + block_i * game->get_N() + block_j)->target / game->get_N(), (game->square_list + block_i * game->get_N() + block_j)->target % game->get_N(), empty);
    }
    void draw_outer_border(Game* game, long long line_radius, bool nice_boundary = true) {
        for (long long i = game->prim_sqr.begin_y; i != game->prim_sqr.begin_y + game->prim_sqr.size; ++i) {
            for (long long j = -line_radius; j < line_radius; ++j) {
                graph_buffer[i * GRAPH_SIZE + (game->prim_sqr.begin_x + j)] = BLACK; /* left edge */
                graph_buffer[i * GRAPH_SIZE + (game->prim_sqr.begin_x + game->prim_sqr.size + j)] = BLACK; /* right edge */
            }
        }
        long long tmp = nice_boundary * line_radius;
        for (long long i = game->prim_sqr.begin_y - tmp; i != game->prim_sqr.begin_y + game->prim_sqr.size + tmp; ++i) {
            for (long long j = -line_radius; j < line_radius; ++j) {
                graph_buffer[(game->prim_sqr.begin_x + j) * GRAPH_SIZE + i] = BLACK; /* top edge */
                graph_buffer[(game->prim_sqr.begin_x + game->prim_sqr.size + j) * GRAPH_SIZE + i] = BLACK; /* bottom edge */
            }
        }
    }
    void fill_game_region(Game* game, const Hue* hue) {
        for (long long i = 0; i != game->get_N() * game->get_N(); ++i)
            fill_block(game, hue, i / game->get_N(), i % game->get_N(), i == game->get_empty());
    }
    bool move_if_asked(Game *game, ExMessage &msg, Hue &hue, int border_line_radius) {
        if (msg.message == WM_LBUTTONUP && custom::is_in_rect(&msg, game->prim_sqr.begin_x, game->prim_sqr.begin_y, game->prim_sqr.size, game->prim_sqr.size)) {
            long long new_i = (msg.y - game->prim_sqr.begin_y) / (game->prim_sqr.size / game->get_N()), new_j = (msg.x - game->prim_sqr.begin_x) / (game->prim_sqr.size / game->get_N()), delta_i = new_i - game->get_empty() / game->get_N(), delta_j = new_j - game->get_empty() % game->get_N();
            if (delta_i * delta_i + delta_j * delta_j == 1) {
                square::swap(game->square_list + game->get_empty(), game->square_list + game->get_empty() + delta_i * game->get_N() + delta_j);
                square::fill_block(game, &hue, game->get_empty() / game->get_N(), game->get_empty() % game->get_N(), false);
                game->move(delta_i, delta_j);
                mciSendStringW(L"play Files\\Sounds\\knock.wav", NULL, 0, NULL) /* playing asynchronously and without stopping the BGM */;
                square::fill_block(game, &hue, game->get_empty() / game->get_N(), game->get_empty() % game->get_N(), true);
                square::draw_outer_border(game, border_line_radius, false /* not necessary */);
                FlushBatchDraw();
                game->increment_cnt();
                return true;
            }
        }
        return false;
    }
}

namespace triangle {
    std::pair<long double, long double> rotate_point(const std::pair<long double, long double> P, long double angle /* clockwise */) {
        return std::pair<long double, long double>(P.first * cos(angle) + P.second * sin(angle), -P.first * sin(angle) + P.second * cos(angle));
    }
    std::pair<long double, long double> vector_add(const std::pair<long double, long double> A, const std::pair<long double, long double> B) {
        return std::pair<long double, long double>(A.first + B.first, A.second + B.second);
    }
    std::pair<long double, long double> vector_minus(const std::pair<long double, long double> A, const std::pair<long double, long double> B) {
        return std::pair<long double, long double>(A.first - B.first, A.second - B.second);
    }
    bool ray_in_segment(const std::pair<long long, long long> P, const std::pair<long long, long long> A, const std::pair<long long, long long> B) {
        if (A.second == B.second)
            return false;
        if (B.second > A.second)
            return (2 * P.second + 1 - 2 * A.second) * (B.first - A.first) > (2 * P.first + 1 - 2 * A.first) * (B.second - A.second) && P.second >= A.second && P.second < B.second;
        else
            return (2 * P.second + 1 - 2 * A.second) * (B.first - A.first) < (2 * P.first + 1 - 2 * A.first) * (B.second - A.second) && P.second < A.second&& P.second >= B.second;
    }
    const long long ASSIGN_PIXEL_STATUS = 2, LEAVE_EMPTY = 1; /* similar to a macro */
    //void fill_block_target; //TODO Maybe imitate that of square::
    void fill_block(Game* game, const Hue* hue, long long* pixel_status, long long N, long long layer, long long k, short mode /* see the line above for macros */) {
        /** Make a copy of vertices and sort them in custom order to find min / max */
        std::pair<long long, long long>* vertices = new std::pair<long long, long long>[3];
        for (long long i = 0; i != 3; ++i)
            vertices[i] = game->triangle_list[layer][k].vertices[i];
        std::sort(vertices, vertices + 3, custom::comp_x);
        long long j_begin = vertices[0].first, j_end = vertices[2].first;
        std::sort(vertices, vertices + 3, custom::comp_y);
        /** ----------------------------------------------------------------------- */
        for (long long i = vertices[0].second; i != vertices[2].second; ++i) {
            for (long long j = j_begin; j != j_end; ++j) {
                long long cnt = 0;
                for (short k = 0; k != 3; ++k)
                    cnt += ray_in_segment(std::pair<long long, long long>(j, i), vertices[k], vertices[(k + 1) % 3]);
                if (cnt == 1) {
                    /** The center of pixel (i, j) lies in the interior of the triangle */
                    if (mode & ASSIGN_PIXEL_STATUS)
                        pixel_status[i * GRAPH_SIZE + j] = k;
                    else if (mode & LEAVE_EMPTY)
                        graph_buffer[i * GRAPH_SIZE + j] = WHITE;
                    else {
                        long long target_k = game->triangle_list[layer][k].target;
                        long long tmp_dir = game->triangle_list[layer][k].direction;
                        if (game->triangle_list[layer][k].type != game->triangle_list[layer][target_k].type)
                            game->triangle_list[layer][k].vertices[tmp_dir] = std::pair<long long, long long>(game->triangle_list[layer][k].vertices[(tmp_dir + 1) % 3].first + game->triangle_list[layer][k].vertices[(tmp_dir + 2) % 3].first - game->triangle_list[layer][k].vertices[tmp_dir].first, game->triangle_list[layer][k].vertices[(tmp_dir + 1) % 3].second + game->triangle_list[layer][k].vertices[(tmp_dir + 2) % 3].second - game->triangle_list[layer][k].vertices[tmp_dir].second);
                        std::pair<long double, long double> k_centroid((game->triangle_list[layer][k].vertices[0].first + game->triangle_list[layer][k].vertices[1].first + game->triangle_list[layer][k].vertices[2].first) / 3.0L, (game->triangle_list[layer][k].vertices[0].second + game->triangle_list[layer][k].vertices[1].second + game->triangle_list[layer][k].vertices[2].second) / 3.0F); /* reflected if necessary */
                        std::pair<long double, long double> target_k_centroid((game->triangle_list[layer][target_k].vertices[0].first + game->triangle_list[layer][target_k].vertices[1].first + game->triangle_list[layer][target_k].vertices[2].first) / 3.0L, (game->triangle_list[layer][target_k].vertices[0].second + game->triangle_list[layer][target_k].vertices[1].second + game->triangle_list[layer][target_k].vertices[2].second) / 3.0F);
                        long double rotate_angle = 3.1415926535897932384626433832795L * 4 * tmp_dir / 3;
                        std::pair<long double, long double> point(j + 0.5L, i + 0.5L);
                        if (game->triangle_list[layer][k].type != game->triangle_list[layer][target_k].type)
                            point = std::pair<long double, long double>(game->triangle_list[layer][k].vertices[(tmp_dir + 1) % 3].first + game->triangle_list[layer][k].vertices[(tmp_dir + 2) % 3].first - point.first, game->triangle_list[layer][k].vertices[(tmp_dir + 1) % 3].second + game->triangle_list[layer][k].vertices[(tmp_dir + 2) % 3].second - point.second);
                        std::pair<long long, long long> target_point = vector_minus(vector_add(rotate_point(vector_minus(point, k_centroid), rotate_angle), target_k_centroid), std::pair<long double, long double>(0.5L, 0.5L));
                        if (game->triangle_list[layer][k].type != game->triangle_list[layer][target_k].type)
                            game->triangle_list[layer][k].vertices[tmp_dir] = std::pair<long long, long long>(game->triangle_list[layer][k].vertices[(tmp_dir + 1) % 3].first + game->triangle_list[layer][k].vertices[(tmp_dir + 2) % 3].first - game->triangle_list[layer][k].vertices[tmp_dir].first, game->triangle_list[layer][k].vertices[(tmp_dir + 1) % 3].second + game->triangle_list[layer][k].vertices[(tmp_dir + 2) % 3].second - game->triangle_list[layer][k].vertices[tmp_dir].second);
                        graph_buffer[i * GRAPH_SIZE + j] = color_hue_BGR(target_point.second, target_point.first, hue);
                    }
                    /** --------------------------------------------------------------- */
                }
            }
        }
    }
    const long long ALGO_DIST_TO_SIDE = 0, ALGO_SMALLER_AND_BIGGER_INTERIOR = 1; /* similar to a macro */
    void draw_outer_border(Game* game, long long N, long long line_radius /* distance is calculated to the center of the pixel instead of the point in the pixel farthest away from the line */, short algo = ALGO_SMALLER_AND_BIGGER_INTERIOR /* see the line above for macros */) {
        std::pair<long long, long long> vertices[3];
        for (long long i = 0; i != 3; ++i)
            vertices[i] = game->prim_tri.vertices[i];
        if (algo == ALGO_SMALLER_AND_BIGGER_INTERIOR) {
            std::pair<long double, long double> vertices_pm[2][3];
            long long j_begin[2], i_begin[2], j_end[2], i_end[2];
            for (short k = 0; k != 3; ++k) {
                long double x = vertices[k].first - (vertices[0].first + vertices[1].first + vertices[2].first) / 3.0L, y = vertices[k].second - (vertices[0].second + vertices[1].second + vertices[2].second) / 3.0L;
                std::pair<long double, long double> delta(line_radius * x * sqrtl(3.0L / (x * x + y * y)), line_radius * y * sqrtl(3.0L / (x * x + y * y)));
                vertices_pm[0][k] = vector_add(vertices[k], delta), vertices_pm[1][k] = vector_minus(vertices[k], delta);
            }
            for (short k = 0; k != 2; ++k)
                j_begin[k] = min_element(vertices_pm[k], vertices_pm[k] + 3, custom::comp_x)->first, j_end[k] = max_element(vertices_pm[k], vertices_pm[k] + 3, custom::comp_x)->first, i_begin[k] = min_element(vertices_pm[k], vertices_pm[k] + 3, custom::comp_y)->second, i_end[k] = max_element(vertices_pm[k], vertices_pm[k] + 3, custom::comp_y)->second;
            for (long long i = i_begin[0]; i != i_end[0]; ++i)
                for (long long j = j_begin[0]; j != j_end[0]; ++j) {
                    long long cnt[2];
                    for (short k = 0; k != 2; ++k) {
                        cnt[k] = 0;
                        if (i >= i_begin[k] && i < i_end[k] && j >= j_begin[k] && j < j_end[k]) {
                            cnt[k] += ray_in_segment(std::pair<long long, long long>(j, i), vertices_pm[k][0], vertices_pm[k][1]);
                            cnt[k] += ray_in_segment(std::pair<long long, long long>(j, i), vertices_pm[k][1], vertices_pm[k][2]);
                            cnt[k] += ray_in_segment(std::pair<long long, long long>(j, i), vertices_pm[k][2], vertices_pm[k][0]);
                        }
                    }
                    if (cnt[0] == 1 && cnt[1] != 1 /* the center of pixel (i, j) lies in the interior of the enlarged triangle but not in the interior of the shirnked triangle */)
                        graph_buffer[i * GRAPH_SIZE + j] = BLACK;
                }
        }
        if (algo == ALGO_DIST_TO_SIDE) {
            for (long long i = 0; i != GRAPH_SIZE; ++i)
                for (long long j = 0; j != GRAPH_SIZE; ++j) {
                    long double dist_sidelines[3];
                    short K = 0, k = 0;
                    do {
                        long long Ax = 2 * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * vertices[(k + 1) % 3].first + (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (2 * j + 1) + (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (2 * vertices[(k + 1) % 3].second - 2 * i - 1), Ay = (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (2 * i + 1) + 2 * (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * vertices[(k + 1) % 3].second + (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (2 * vertices[(k + 1) % 3].first - 2 * j - 1), Bx = 2 * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * vertices[(k + 2) % 3].first + (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (2 * j + 1) + (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (2 * vertices[(k + 2) % 3].second - 2 * i - 1), By = (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (2 * i + 1) + 2 * (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * vertices[(k + 2) % 3].second + (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (2 * vertices[(k + 2) % 3].first - 2 * j - 1), delta = (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) + (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second);
                        if ((delta * (2 * i + 1) - Ay) * (delta * (2 * i + 1) - By) <= 0 && (delta * (2 * j + 1) - Ax) * (delta * (2 * j + 1) - Bx) <= 0 /* the center of the pixel is close to the side of the triangle besides being close to the extended sideline */)
                            dist_sidelines[K++] = abs((vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (j + 0.5L) - (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (i + 0.5L) + (vertices[(k + 2) % 3].first * vertices[(k + 1) % 3].second - vertices[(k + 1) % 3].first * vertices[(k + 2) % 3].second)) / sqrtl((vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) + (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second)) /* make sure it is floating point division */;
                    } while (++k != 3);
                    if (K && *std::min_element(dist_sidelines, dist_sidelines + K) < line_radius)
                        graph_buffer[i * GRAPH_SIZE + j] = BLACK;
                }
        }

    }
    void fill_game_region(Game* game, const Hue* hue, long long *pixel_status, bool assign_pixel_status, long long empty = -1 /* valid only if assign_pixel_status is false */) {
        for (long long i = 0; i != custom::raise_to_power(4, game->get_N() - 1); ++i)
            fill_block(game, hue, pixel_status, game->get_N(), game->get_N() - 1, i, assign_pixel_status? ASSIGN_PIXEL_STATUS: (i == empty? LEAVE_EMPTY: 0));
    }
}


void play_bgm(short bgm /* 0 or 1 */) {
    PlaySoundW((L"Files\\Music\\" + std::to_wstring(bgm) + L".wav").c_str(), NULL, SND_ASYNC | SND_LOOP | SND_FILENAME); //TODO need to avoid absolute path but this function needs absolute path? see Lv Xinyao's project
}

void play_rand_bgm() {
    play_bgm(current_bgm = custom::rand(2));
}


/** ----- ARCHIVED -----
*   void write_err_log(std::string err) {
*       std::ofstream err_log;
*       err_log.open("Files\\ERR_LOG.txt", std::fstream::app);
*       err_log << err;
*       err_log.close();
*   }
*/

void fill_entire_graph(const Hue* hue) {
    /** Fill the entire graph by a hue of color */
    for (long long i = 0; i != GRAPH_SIZE; ++i) {
        for (long long j = 0; j != GRAPH_SIZE; ++j) {
            graph_buffer[i * GRAPH_SIZE + j] = color_hue_BGR(i, j, hue);
        }
    }
    /** --------------------------------------- */
}

//void game(short mode, long long N /* positive integer, made greater than 1 for something to play */, const Hue* hue) {
//    if (mode == SQUARE_MODE) {
//
//        //FlushBatchDraw();//TODO Learn where it should be put at
//        /** Awaiting for mouse action while not winned */
//        ExMessage msg;
//        while (!sqrgame::check_solved(square_list, N)) {
//            peekmessage(&msg, EX_MOUSE); //TODO peekmessage returns a boolean value, need to use or not?
//            if (msg.message == WM_LBUTTONUP) {
//                if (msg.x >= prim_sqr.begin_x && msg.x < prim_sqr.begin_x + prim_sqr.size && msg.y >= prim_sqr.begin_y && msg.y < prim_sqr.begin_y + prim_sqr.size) {
//                    long long new_i = (msg.y - prim_sqr.begin_y) / block_size, new_j = (msg.x - prim_sqr.begin_x) / block_size, delta_i = new_i - empty_i, delta_j = new_j - empty_j;
//                    if (delta_i * delta_i + delta_j * delta_j == 1) {
//                        sqrgame::swap(square_list + empty_i * N + empty_j, square_list + new_i * N + new_j);
//                        sqrgame::fill_block(graph_buffer, square_list, hue, N, empty_i, empty_j, false);
//                        empty_i += delta_i, empty_j += delta_j;
//                        //TODO Sound effects
//                        sqrgame::fill_block(graph_buffer, square_list, hue, N, empty_i, empty_j, true);
//                        sqrgame::draw_outer_border(graph_buffer, N, 3, false /* not necessary */);
//                    }
//                }
//            }
//        }
//        /** ------------------------------------------ */
//    }
//    else if (mode == TRIANGLE_MODE) {
//        long long* pixel_status = new long long [GRAPH_SIZE * GRAPH_SIZE];
//        for (long long i = 0; i != GRAPH_SIZE * GRAPH_SIZE; ++i)
//            pixel_status[i] = -1;
//        trigame::triangle_construction(graph_buffer, triangle_list, N - 1, 0, 0);
//        trigame::fill_game_region(graph_buffer, triangle_list, hue, pixel_status, N, true);
//        long long empty = custom::rand(custom::raise_to_power(4, N - 1));
//        /** Shuffle the game region */
//        for (long long k = 0; k != shuffle_times; ++k) { //TODO Need to consider the case where no shuffle is done
//            long long tmp_rand = custom::rand(3);
//            std::pair<long double, long double> empty_centroid((triangle_list[N - 1][empty].vertices[0].first + triangle_list[N - 1][empty].vertices[1].first + triangle_list[N - 1][empty].vertices[2].first) / 3.0L, (triangle_list[N - 1][empty].vertices[0].second + triangle_list[N - 1][empty].vertices[1].second + triangle_list[N - 1][empty].vertices[2].second) / 3.0F);
//            std::pair<long long, long long> target_centroid(2.0L * empty_centroid.first - triangle_list[N - 1][empty].vertices[tmp_rand].first - 0.5L, 2.0L * empty_centroid.second - triangle_list[N - 1][empty].vertices[tmp_rand].second - 0.5L);
//            //TODO the well behavior of this routine depends somehow on the precision
//            long long target;
//            if (target_centroid.first >= 0 && target_centroid.first < GRAPH_SIZE && target_centroid.second >= 0 && target_centroid.second < GRAPH_SIZE && (target = pixel_status[target_centroid.second * GRAPH_SIZE + target_centroid.first]) != -1) {
//                if (target == empty /* I believe this should not happen */)
//                    throw std::runtime_error("Should not equal. empty = " + std::to_string(empty) + "; target = " + std::to_string(target) + ".\n");
//                else if (!trigame::check_adjacent_and_swap(triangle_list, N, target, empty)) /* I believe the line below should not execute */ {
//                    throw std::runtime_error("Not adjacent. empty = " + std::to_string(empty) + "; target = " + std::to_string(target));
//                }
//                empty = target;
//            }
//        }
//        /** ----------------------- */
//        trigame::fill_game_region(graph_buffer, triangle_list, hue, pixel_status, N, false, empty);
//        trigame::draw_outer_border(graph_buffer, N, 2);
//        /** Awaiting for mouse action while not winned */
//        ExMessage msg;
//        while (!trigame::check_solved(triangle_list, N)) {
//            peekmessage(&msg, EX_MOUSE);
//            if (msg.message == WM_LBUTTONUP) {
//                long long target;
//                if ((target = pixel_status[msg.y * GRAPH_SIZE + msg.x]) != -1 && target != empty && trigame::check_adjacent_and_swap(triangle_list, N, target, empty)) {
//                    trigame::fill_block(graph_buffer, triangle_list, hue, pixel_status, N, N - 1, empty, 0);
//                    bool need_draw_border = triangle_list[N - 1][empty].boundary > 0;
//                    empty = target;
//                    //TODO Sound effects
//                    trigame::fill_block(graph_buffer, triangle_list, hue, pixel_status, N, N - 1, empty, trigame::LEAVE_EMPTY);
//                    if (need_draw_border || triangle_list[N - 1][empty].boundary > 0)
//                        trigame::draw_outer_border(graph_buffer, N, 2);
//                }
//            }
//        }
//        /** ------------------------------------------ */
//    }
//    /** After winned */
//    PlaySoundW(NULL, 0, 0); /* stop bgm */
//    std::this_thread::sleep_for(std::chrono::milliseconds(200));
//    PlaySoundW(L"Files\\Sounds\\spalsh.wav", NULL, SND_ASYNC | SND_FILENAME); /* about 350 milliseconds delay */
//    std::this_thread::sleep_for(std::chrono::milliseconds(30));
//    fill_entire_graph(graph_buffer, hue);
//    std::this_thread::sleep_for(std::chrono::milliseconds(420));
//    play_rand_bgm(); /* play bgm */
//    /** ------------ */
//}

void countdown_for(Game* game, short countdown /* a single digit odd number */) /* this routine includes sleep_for */ { //TODO Implement for Triangle
    long long empty_block_begin_x = game->prim_sqr.begin_x + (game->get_empty(true) % game->get_N()) * (game->prim_sqr.size / game->get_N()), empty_block_begin_y = game->prim_sqr.begin_y + (game->get_empty(true) / game->get_N()) * (game->prim_sqr.size / game->get_N());
    settextcolor(BLACK);
    setbkcolor(WHITE);
    setbkmode(OPAQUE);
    if (game->get_N() <= 4) {
        settextstyle(70, 0, L"仿宋", 0, 0, 300, false, false, false);
        outtextxy(empty_block_begin_x + game->prim_sqr.size / game->get_N() / 2 - 85, empty_block_begin_y + game->prim_sqr.size / game->get_N() / 2 - 85, L"原");
        outtextxy(empty_block_begin_x + game->prim_sqr.size / game->get_N() / 2 + 15, empty_block_begin_y + game->prim_sqr.size / game->get_N() / 2 + 15, L"图");
        settextcolor(RED);
        short bottom_left_top_right = custom::rand(2);
        for (short k = 0; k != countdown; ++k) {
            settextstyle(100, 0, L"Consolas", 0, 0, 400, false, false, false);
            outtextxy(empty_block_begin_x + game->prim_sqr.size / game->get_N() / 2 + (bottom_left_top_right ? 30 : -70), empty_block_begin_y + game->prim_sqr.size / game->get_N() / 2 + (bottom_left_top_right ? -100 : 0), std::to_wstring(countdown - k).c_str());
            settextstyle(100, 0, L"Calibri", 0, 0, 400, false, false, false);
            outtextxy(empty_block_begin_x + game->prim_sqr.size / game->get_N() / 2 + (bottom_left_top_right ? -70 : 10), empty_block_begin_y + game->prim_sqr.size / game->get_N() / 2 + (bottom_left_top_right ? 0 : -100), k % 2 ? L"    " : (bottom_left_top_right ? L"↗" : L"↙"));
            FlushBatchDraw();
            std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        }
    }
    else if (game->get_N() == 5) {
        settextstyle(60, 0, L"仿宋", 0, 0, 400, false, false, false);
        outtextxy(empty_block_begin_x + game->prim_sqr.size / game->get_N() / 2 - 75, empty_block_begin_y + game->prim_sqr.size / game->get_N() / 2 - 75, L"原");
        outtextxy(empty_block_begin_x + game->prim_sqr.size / game->get_N() / 2 + 15, empty_block_begin_y + game->prim_sqr.size / game->get_N() / 2 + 15, L"图");
        settextcolor(RED);
        for (short k = 0; k != countdown; ++k) {
            settextstyle(50, 0, L"Consolas", 0, 0, 400, false, false, false);
            outtextxy(empty_block_begin_x + game->prim_sqr.size / game->get_N() / 2 - 15, empty_block_begin_y + game->prim_sqr.size / game->get_N() / 2 - 20, std::to_wstring(countdown - k).c_str());
            FlushBatchDraw();
            std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        }
    }
    else if (game->get_N() >= 6) {
        settextcolor(RED);
        for (short k = 0; k != countdown; ++k) {
            settextstyle(30, 0, L"Consolas", 0, 0, 400, false, false, false);
            outtextxy(empty_block_begin_x + game->prim_sqr.size / game->get_N() / 2 - custom::consolas_width(30) / 2, empty_block_begin_y + game->prim_sqr.size / game->get_N() / 2 - 15, std::to_wstring(countdown - k).c_str());
            FlushBatchDraw();
            std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        }
    }
    setbkmode(TRANSPARENT);
    flushmessage();
}

void show_timer(const short digits_before_decimal /* positive integer */, const long long milliseconds, bool made_red_and_bold) {
    setlinestyle(PS_DASH, 3);
    setfillcolor(WHITE);
    setlinecolor(0xe2e2e2);
    fillrectangle(712, 15, 920, 65);
    std::wstring out_str;
    long long x = milliseconds / 10000;
    for (short k = 0; k != digits_before_decimal - 1; ++k) {
        if (x)
            out_str = (wchar_t)(L'0' + x % 10) + out_str, x /= 10;
        else
            out_str = L' ' + out_str;
    }
    out_str += (wchar_t)(L'0' + milliseconds / 1000 % 10);
    out_str += L".";
    out_str += (wchar_t)(L'0' + milliseconds / 100 % 10);
    settextcolor(made_red_and_bold ? RED: 0x404040);
    settextstyle(30, 0, L"楷体", 0, 0, made_red_and_bold? 900: 500, false, false, false);
    outtextxy(718, 23, L"时间");
    settextstyle(30, 0, L"Consolas", 0, 0, made_red_and_bold? 900: 500, false, false, false);
    outtextxy(782, 23, (L": " + out_str + L" s").c_str());
    FlushBatchDraw();
}

void sound_after_game(bool failed_winned) {
    PlaySoundW(NULL, 0, 0); /* stop current bgm */
    std::this_thread::sleep_for(std::chrono::milliseconds(failed_winned? 200: 150));
    PlaySoundW(failed_winned? L"Files\\Sounds\\splash.wav": L"Files\\Sounds\\gameover.wav", NULL, SND_ASYNC | SND_FILENAME); /* about 350 or 2250 milliseconds delay */
    std::this_thread::sleep_for(std::chrono::milliseconds(failed_winned? 400: 2250));
    play_rand_bgm(); /* play bgm */
    flushmessage();
}

namespace image_res {
    IMAGE background[2], exit, pause, play, play_small, restart, sound[2], star_empty, star_filled[2];
    void load() {
        for (short i = 0; i != 2; ++i)
            loadimage(&background[i], (L"Files\\Images\\Elements\\Flower_" + std::to_wstring(i) + L".png").c_str(), GRAPH_SIZE / 5 * 2, GRAPH_SIZE / 5 * 2);
        loadimage(&exit, L"Files\\Images\\Elements\\Exit.png", 80, 80);
        loadimage(&pause, L"Files\\Images\\Elements\\Pause.png", 50, 50);
        loadimage(&play, L"Files\\Images\\Elements\\Play.png", 50, 50);
        loadimage(&play_small, L"Files\\Images\\Elements\\Play.png", 24, 24);
        loadimage(&restart, L"Files\\Images\\Elements\\Restart.png", 35, 35);
        for (short i = 0; i != 2; ++i)
            loadimage(&sound[i], i ? L"Files\\Images\\Elements\\Sound_on.png" : L"Files\\Images\\Elements\\Sound_off.png", 50, 50);
        loadimage(&star_empty, L"Files\\Images\\Elements\\Star_empty.png", 25, 25);
        for (short i = 0; i != 2; ++i)
            loadimage(&star_filled[i], (L"Files\\Images\\Elements\\Star_filled" + std::to_wstring(i) + L".png").c_str(), 25, 25);
    }
}

int main() {
    GetCurrentDirectoryW(MAX_PATH, working_path) /* initialization */;
#ifdef DEBUG
    std::cout << "Debugging now.\n";
    initgraph(GRAPH_SIZE, GRAPH_SIZE, EX_SHOWCONSOLE | EX_NOCLOSE | EX_NOMINIMIZE);
#else
    initgraph(GRAPH_SIZE, GRAPH_SIZE, EX_NOCLOSE | EX_NOMINIMIZE);
#endif
    play_rand_bgm();
    setbkmode(TRANSPARENT);
    image_res::load();
    Hue hue_play = { { 0, 0, 0 }, { 0, 256, 0 }, { 256, 0, 0 } };
    Hue hue_sky = { { 0, 0, 255 }, { 0, 0, 0 }, { 0, 0, 0 } } /* some values are to be assigned */;
    graph_buffer = GetImageBuffer();
    ExMessage msg;
    BeginBatchDraw();
    /** Game */ {
        Game* this_game = nullptr /* an initialization is required, otherwise results in fatal error LNK1257 */;
        page = START_PAGE;
        while (true) {
            flushmessage();
            /* it is necessary to deal with START_PAGE and EXIT_PAGE cases first  */
            if (page == START_PAGE) {
#ifdef DEBUG
                std::cout << "page is START_PAGE now.\n";
#endif
                Hue hue_play_small = { { 64, 64, 0 }, { 0, 1280, 0 }, { 1280, 0, 0 } } /* the size of the small square is GRAPH_SIZE / 10 */;
                Game small_square_game(SQUARE_MODE, 2);
                small_square_game.set_prim_sqr(GRAPH_SIZE / 5, GRAPH_SIZE / 5, GRAPH_SIZE / 10);
                small_square_game.set_sqr(true);
                /** Shuffle the small game so that one move is required to solve it */
                do {
                    small_square_game.shuffle(1);
                } while (small_square_game.check_solved());
                /** --------------------------------------------------------------- */
                small_square_game.reset_timer() /* maybe not used */;
                small_square_game.reset_cnt() /* maybe not used */;
                const long long countdown_before_exit = 3 /* a single digit number */;
                long long time_elapsed_since_exit_request /* in seconds */;
                std::chrono::milliseconds flash_period(current_bgm ? 255 : 170 /* a BGM specific value figured out by rough testing */), hint_appear(600 + custom::rand(600));
                std::chrono::steady_clock::time_point begin_time = std::chrono::steady_clock::now(), last_time = begin_time, exit_request_time;
                bool hue_sky_up_down = false, bg_refresh = true, hint_enabled = false, play_icon_enabled = true /* if hint is enabled */, refresh_after_winned = false /* refresh for the last time */, prepare_exit = false;
                int refresh_cycle_cnt = 0;
                while (true) {
                    std::chrono::steady_clock::time_point time_this_loop = std::chrono::steady_clock::now();
                    if (prepare_exit) {
                        time_elapsed_since_exit_request = std::chrono::duration_cast<std::chrono::seconds>(time_this_loop - exit_request_time).count();
                        if (time_elapsed_since_exit_request >= countdown_before_exit) {
                            page = EXIT_PAGE;
                            break;
                        }
                    }
                    hint_enabled = hint_enabled || std::chrono::duration_cast<std::chrono::milliseconds>(time_this_loop - begin_time) >= hint_appear;
                    if (bg_refresh) {
                        /** One refresh cycle */
                        hue_sky.top_left_RGB[0] = hue_sky_up_down ? 127 : 80;
                        hue_sky.top_left_RGB[1] = hue_sky_up_down ? 223 : 192;
                        hue_sky.total_change_i[0] = hue_sky_up_down ? -48 : 48;
                        hue_sky.total_change_i[1] = hue_sky_up_down ? -32 : 32;
                        for (int i = 0; i != GRAPH_SIZE; ++i)
                            for (int j = 0; j != GRAPH_SIZE; ++j)
                                graph_buffer[i * GRAPH_SIZE + j] = color_hue_BGR(i, j, &hue_sky);
                        for (short i = 0; i != 2; ++i)
                            putimage(GRAPH_SIZE / 5 * 3, 0, &image_res::background[i], i ? SRCPAINT : SRCAND) /* successful but picture jagged */;
                        putimage(5, 0, &image_res::exit, SRCAND);
                        setlinestyle(PS_DASH, 20);
                        setlinecolor(BLACK);
                        circle(500, 500, 180);
                        if (prepare_exit) {
                            settextcolor(RED);
                            settextstyle(150, 0, L"Consolas", 0, 0, 400, false, false, false);
                            outtextxy(500 - custom::consolas_width(150) / 2, 425, std::to_wstring(countdown_before_exit - time_elapsed_since_exit_request /* a single digit */).c_str());
                            //TODO exit.png下面弄退出的文字提示; 或者3 2 1的底下“Exit?”
                            FlushBatchDraw();
                        }
                        else {
                            settextstyle(100, 0, L"仿宋", 0, 0, 300, false, false, false);
                            settextcolor(WHITE);
                            outtextxy(385, 385, L"滑");
                            outtextxy(515, 515, L"动");
                            outtextxy(385, 515, L"拼");
                            outtextxy(515, 385, L"图");
                            setlinestyle(PS_DASH, 5);
                            line(385, 500, 614, 500);
                            line(500, 385, 500, 614);
                            if (hint_enabled) {
                                settextcolor(0xe2e2 /* gray */);
                                if (refresh_cycle_cnt % 3 != 0 && play_icon_enabled) {
                                    settextstyle(30, 0, L"Calibri", 0, 0, 600, false, true, false);
                                    outtextxy(207, 160, L"Click it!");
                                }
                                if (refresh_cycle_cnt % 9 != 0) {
                                    settextstyle(20, 0, L"Candara", 0, 0, 400, false, false, false);
                                    outtextxy(204, 308, L"(sliding block)");
                                }
                            }
                            square::fill_game_region(&small_square_game, &hue_play_small);
                            if (refresh_after_winned) {
                                square::fill_block(&small_square_game, &hue_play_small, small_square_game.get_empty() / small_square_game.get_N(), small_square_game.get_empty() % small_square_game.get_N(), false);
                            }
                            square::draw_outer_border(&small_square_game, 2, true);
                            if (refresh_after_winned) {
                                FlushBatchDraw();
                                sound_after_game(true);
                                page = SQUARE_MODE | HUE_MODE | GAME_PAGE;
                                break;
                            }
                            else {
                                if (hint_enabled && play_icon_enabled && refresh_cycle_cnt % 3 != 0) {
                                    int target = small_square_game.square_list[small_square_game.get_empty()].target;
                                    putimage(small_square_game.prim_sqr.begin_x + (target % 2) * (small_square_game.prim_sqr.size / 2) + small_square_game.prim_sqr.size / 4 - 12, small_square_game.prim_sqr.begin_y + (target / 2) * (small_square_game.prim_sqr.size / 2) + small_square_game.prim_sqr.size / 4 - 12, &image_res::play_small, SRCAND) /* the original idea is to draw a circle with color 0xe2e2e2 */;
                                }
                                FlushBatchDraw();
                            }
                        }
                        /* ------------------ */
                        hue_sky_up_down = !hue_sky_up_down, bg_refresh = false, refresh_cycle_cnt++;
                    }
                    else if (!bg_refresh) {
                        if (std::chrono::duration_cast<std::chrono::milliseconds>(time_this_loop - last_time) >= flash_period)
                            last_time = time_this_loop, bg_refresh = true;
                    }
                    bool peeked = peekmessage(&msg, EX_MOUSE);
                    if (peeked) {
                        if (!prepare_exit && square::move_if_asked(&small_square_game, msg, hue_play_small, 2)) {
                            if (small_square_game.check_solved()) {
                                refresh_after_winned = true /* after small game solved */, bg_refresh = true /* refresh immediately within the next loop */;
                                continue;
                            }
                            else {
                                long long empty = small_square_game.get_empty(), target = small_square_game.square_list[empty].target;
                                bool only_one_step_till_solved = true;
                                for (long long i = 0; only_one_step_till_solved && i != small_square_game.get_N() * small_square_game.get_N(); ++i) {
                                    if (i != empty && i != target)
                                        only_one_step_till_solved = only_one_step_till_solved && small_square_game.square_list[i].target == i;
                                }
                                play_icon_enabled = only_one_step_till_solved;
                            }
                        }
                        if (!prepare_exit && custom::is_in_rect(&msg, 23, 24, 27, 32)) {
                            prepare_exit = true, exit_request_time = time_this_loop, bg_refresh = true /* refresh immediately within the next loop */;
                            continue;
                        }
                        if (prepare_exit && !custom::is_in_rect(&msg, 23, 24, 27, 32)) {
                            prepare_exit = false, bg_refresh = true /* refresh immediately within the next loop */;
                            continue;
                        }
                        if (msg.message == WM_LBUTTONUP && custom::is_in_rect(&msg, 23, 24, 27, 32)) {
                            page = EXIT_PAGE;
                            break;
                        }
                    }
                }
                continue;
            }
            else if (page == EXIT_PAGE) {
#ifdef DEBUG
                std::cout << "page is EXIT_PAGE now.\n";
#endif
                Hue hue_play_small = { { 64, 64, 0 }, { 0, 1280, 0 }, { 1280, 0, 0 } } /* the size of the small square is GRAPH_SIZE / 10 */;
                std::chrono::milliseconds flash_period(current_bgm ? 255 : 170 /* a BGM specific value figured out by rough testing */), till_exit(1000);
                std::chrono::steady_clock::time_point begin_time = std::chrono::steady_clock::now(), last_time = begin_time;
                bool hue_sky_up_down = false, bg_refresh = true;
                while (true) {
                    std::chrono::steady_clock::time_point time_this_loop = std::chrono::steady_clock::now();
                    if (std::chrono::duration_cast<std::chrono::milliseconds>(time_this_loop - begin_time) >= till_exit) {
                        break;
                    }
                    if (bg_refresh) {
                        /** One refresh cycle */
                        hue_sky.top_left_RGB[0] = hue_sky_up_down ? 127 : 80;
                        hue_sky.top_left_RGB[1] = hue_sky_up_down ? 223 : 192;
                        hue_sky.total_change_i[0] = hue_sky_up_down ? -48 : 48;
                        hue_sky.total_change_i[1] = hue_sky_up_down ? -32 : 32;
                        for (int i = 0; i != GRAPH_SIZE; ++i)
                            for (int j = 0; j != GRAPH_SIZE; ++j)
                                graph_buffer[i * GRAPH_SIZE + j] = color_hue_BGR(i, j, &hue_sky);
                        for (short i = 0; i != 2; ++i)
                            putimage(GRAPH_SIZE / 5 * 3, 0, &image_res::background[i], i ? SRCPAINT : SRCAND) /* successful but picture jagged */;
                        setlinestyle(PS_DASH, 20);
                        setlinecolor(BLACK);
                        circle(500, 500, 180);
                        settextstyle(100, 0, L"仿宋", 0, 0, 300, false, false, false);
                        settextcolor(WHITE);
                        outtextxy(385, 385, L"再");
                        outtextxy(515, 515, L"见");
                        setlinestyle(PS_DASH, 5);
                        line(385, 500, 614, 500);
                        line(500, 385, 500, 614);
                        FlushBatchDraw();
                        /* ------------------ */
                        hue_sky_up_down = !hue_sky_up_down, bg_refresh = false;
                    }
                    else if (!bg_refresh) {
                        if (std::chrono::duration_cast<std::chrono::milliseconds>(time_this_loop - last_time) >= flash_period)
                            last_time = time_this_loop, bg_refresh = true;
                    }
                }
                break /* exit */;
            }
            else if ((page & SQR_TRI_FILTER) == SQUARE_MODE && (page & HUE_IMG_FILTER) == HUE_MODE && (page & PAGE_TYPE) == PAUSED_PAGE) {
#ifdef DEBUG
                std::cout << "page is SQUARE_MODE | HUE_MODE | PAUSED_PAGE now.\n";
#endif
                /** Construct a new game for user to test while real game paused */
                Game new_tmp_game(this_game->get_mode(), this_game->get_N());
                new_tmp_game.set_prim_sqr(GRAPH_SIZE / 2, GRAPH_SIZE / 2, 350);
                new_tmp_game.set_sqr(true);
                new_tmp_game.set_original_empty_block(this_game->get_empty(true));
                new_tmp_game.reset_timer() /* maybe not used */;
                new_tmp_game.reset_cnt();
                /* ------------------------------------------------------------- */
                const long long countdown_before_exit = 3 /* a single digit number */, countdown_before_resume = 3 /* a single digit number */;
                long long time_elapsed_since_exit_request /* in seconds */, time_elapsed_since_resume_request /* in seconds */;
                std::chrono::milliseconds flash_period(current_bgm ? 255 : 170 /* a BGM specific value figured out by rough testing */);
                std::chrono::steady_clock::time_point last_time = std::chrono::steady_clock::now(), exit_request_time, resume_request_time;
                bool hue_sky_up_down = false, bg_refresh = true, prepare_exit = false, prepare_resume = false /* they must not be true at the same time */, ignore_mouse_hover_on_resume = true /* if not ignored at the first time, the behavior is strange when PAUSED_PAGE is entered by hovering the mouse on the pause sign in GAME_PAGE */;
                while (true) {
                    std::chrono::steady_clock::time_point time_this_loop = std::chrono::steady_clock::now();
                    if (prepare_exit) {
                        time_elapsed_since_exit_request = std::chrono::duration_cast<std::chrono::seconds>(time_this_loop - exit_request_time).count();
                        if (time_elapsed_since_exit_request >= countdown_before_exit) {
                            page = START_PAGE;
                            break;
                        }
                    }
                    else if (prepare_resume) {
                        time_elapsed_since_resume_request = std::chrono::duration_cast<std::chrono::seconds>(time_this_loop - resume_request_time).count();
                        if (time_elapsed_since_resume_request >= countdown_before_resume) {
                            page = SQUARE_MODE | HUE_MODE | RESUMED | GAME_PAGE;
                            break;
                        }
                    }
                    if (bg_refresh) {
                        hue_sky.top_left_RGB[0] = hue_sky_up_down ? 127 : 80;
                        hue_sky.top_left_RGB[1] = hue_sky_up_down ? 223 : 192;
                        hue_sky.total_change_i[0] = hue_sky_up_down ? -48 : 48;
                        hue_sky.total_change_i[1] = hue_sky_up_down ? -32 : 32;
                        for (int i = 0; i != GRAPH_SIZE; ++i)
                            for (int j = 0; j != GRAPH_SIZE; ++j)
                                graph_buffer[i * GRAPH_SIZE + j] = color_hue_BGR(i, j, &hue_sky);
                        putimage(5, 0, &image_res::exit, SRCAND);
                        putimage(85, 15, &image_res::play, SRCAND); //LEARN Why the image will flash if put after drawing of square and right before flusing
                        this_game->set_prim_sqr(80, 80, prepare_resume? 840: 350);
                        this_game->set_sqr(false);
                        square::fill_game_region(this_game, &hue_play);
                        square::draw_outer_border(this_game, 3, true);
                        if (prepare_resume || prepare_exit) {
                            settextcolor(RED);
                        }
                        if (prepare_resume) {
                            long long empty_block_begin_x = this_game->prim_sqr.begin_x + (this_game->get_empty() % this_game->get_N()) * (this_game->prim_sqr.size / this_game->get_N()), empty_block_begin_y = this_game->prim_sqr.begin_y + (this_game->get_empty() / this_game->get_N()) * (this_game->prim_sqr.size / this_game->get_N());
                            long long text_height;
                            switch (this_game->get_N()) {
                            case 1:
                            case 2:
                            case 3:
                            case 4:
                                text_height = 100;
                                break;
                            case 5:
                                text_height = 50;
                                break;
                            default:
                                text_height = 30;
                            } //TODO More testing
                            setbkmode(OPAQUE);
                            setbkcolor(WHITE);
                            settextstyle(text_height, 0, L"Consolas", 0, 0, 400, false, false, false);
                            long long tmp_x = empty_block_begin_x + this_game->prim_sqr.size / this_game->get_N() / 2, tmp_y = empty_block_begin_y + this_game->prim_sqr.size / this_game->get_N() / 2;
                            outtextxy(tmp_x - custom::consolas_width(text_height) / 2, tmp_y - text_height / 2, std::to_wstring(countdown_before_resume - time_elapsed_since_resume_request /* a single digit */).c_str());
                            settextstyle(20, 0, L"Candara", 0, 0, 400, false, false, false);
                            outtextxy(tmp_x - 43, tmp_y + text_height / 2 + 3, L"- Going back -");
                            square::draw_outer_border(this_game, 3, true);
                            setbkmode(TRANSPARENT);
                        }
                        else if (prepare_exit) {
                            settextstyle(150, 0, L"Consolas", 0, 0, 400, false, false, false);
                            outtextxy(500 - custom::consolas_width(150) / 2, 425, std::to_wstring(countdown_before_exit - time_elapsed_since_exit_request /* a single digit */).c_str());
                            settextstyle(40, 0, L"Candara", 0, 0, 400, false, false, false);
                            outtextxy(500 - 75, 425 + 150 + 3, L"- Quitting -");
                        }
                        else {
                            for (long long i = 0; i != this_game->get_N() * this_game->get_N(); ++i)
                                square::fill_block_target(&new_tmp_game, &hue_play, i / this_game->get_N(), i % this_game->get_N(), i / this_game->get_N(), i % this_game->get_N(), i == new_tmp_game.get_empty());
                            square::draw_outer_border(&new_tmp_game, 3, true);
                            long long tmp;
                            show_timer(3, tmp = std::chrono::duration_cast<std::chrono::milliseconds>(this_game->get_timer()).count(), tmp * 10 > std::chrono::duration_cast<std::chrono::milliseconds>(this_game->get_time_limit()).count() * 9);
                        }
                        FlushBatchDraw();
                        hue_sky_up_down = !hue_sky_up_down, bg_refresh = false;
                    }
                    else {
                        if (std::chrono::duration_cast<std::chrono::milliseconds>(time_this_loop - last_time) >= flash_period)
                            last_time = time_this_loop, bg_refresh = true;
                    }
                    bool peeked = peekmessage(&msg, EX_MOUSE | EX_KEY);
                    if (peeked) {
                        if (msg.message == WM_KEYUP /* suppose the keyup after keydown is successfully detected */ && msg.vkcode == VK_ESCAPE) {
                            page = EXIT_PAGE;
                            break;
                        } //TODO Escape in almost all pages; documentate
                        if (msg.message == WM_MOUSEMOVE /* I believe this suffices */)
                            ignore_mouse_hover_on_resume = ignore_mouse_hover_on_resume && custom::is_in_circ(&msg, 110, 38, 21) /* it is ok to continue the flow using msg possibly again */;
                        if (prepare_resume) {
                            if (msg.message == WM_KEYUP && msg.vkcode == VK_BACK || msg.message == WM_LBUTTONUP && custom::is_in_circ(&msg, 110, 38, 21)) {
                                page = SQUARE_MODE | HUE_MODE | RESUMED | GAME_PAGE;
                                break;
                            }
                            if (!custom::is_in_circ(&msg, 110, 38, 21)) {
                                prepare_resume = false;
                                continue;
                            }
                        }
                        else if (prepare_exit) {
                            if (msg.message == WM_LBUTTONUP && custom::is_in_rect(&msg, 23, 24, 27, 32)) {
                                page = START_PAGE;
                                break;
                            }
                            if (!custom::is_in_rect(&msg, 23, 24, 27, 32)) {
                                prepare_exit = false;
                                continue;
                            }
                        }
                        else {
                            if (!square::move_if_asked(&new_tmp_game, msg, hue_play, 3)) {
                                prepare_exit = custom::is_in_rect(&msg, 23, 24, 27, 32);
                                prepare_resume = !prepare_exit && !ignore_mouse_hover_on_resume && (msg.message == WM_KEYDOWN && msg.vkcode == VK_SPACE || custom::is_in_circ(&msg, 110, 38, 21));
                                if (prepare_exit)
                                    exit_request_time = time_this_loop, bg_refresh = true /* refresh immediately within the next loop */;
                                if (prepare_resume)
                                    resume_request_time = time_this_loop, bg_refresh = true /* refresh immediately within the next loop */;
                            }
                        }
                    }
                }
                continue;
            }
            else if ((page & SQR_TRI_FILTER) == SQUARE_MODE && (page & HUE_IMG_FILTER) == HUE_MODE && (page & PAGE_TYPE) == GAME_PAGE) {
#ifdef DEBUG
                if (page & RESUMED_FILTER)
                    std::cout << "page is SQUARE_MODE | HUE_MODE | GAME_PAGE | RESUMED now.\n";
                else
                    std::cout << "page is SQUARE_MODE | HUE_MODE | GAME_PAGE now.\n";
#endif
                fill_entire_graph(&hue_play);
                if (!(page & RESUMED)) {
                    /** Initialize a game */
                    Game game(SQUARE_MODE, 9) /* this should not be cleared as long as main() is not returned, so a pointer to it is valid */; //TODO Setting change this constant, no larger than 9
                    game.set_prim_sqr(80, 80, 840 /* will be modified by the program to be a multiple of N */);
                    /** Shuffle the small game so that at least one move is required to solve it */
                    do {
                        game.set_sqr(true);
                        game.shuffle(SHUFFLE_TIMES);
                    } while (game.check_solved());
                    /** ------------------------------------------------------------------------ */
                    for (long long i = 0; i != game.get_N() * game.get_N(); ++i)
                        square::fill_block_target(&game, &hue_play, i / game.get_N(), i % game.get_N(), i / game.get_N(), i % game.get_N(), i == game.get_empty(true));
                    square::draw_outer_border(&game, 3, true);
                    countdown_for(&game, 3);
                    game.reset_timer();
                    game.reset_cnt();
                    game.set_time_limit(std::chrono::seconds(100)); //TODO Customize
                    /** ----------------- */
                    this_game = &game;
                }
                square::fill_game_region(this_game, &hue_play);
                square::draw_outer_border(this_game, 3, (page & RESUMED_FILTER) == RESUMED);
                putimage(5, 0, &image_res::exit, MERGEPAINT);
                putimage(85, 15, &image_res::pause, MERGEPAINT);
                FlushBatchDraw();
                const long long countdown_before_exit = 3 /* a single digit number */, countdown_before_pause = 3 /* a single digit number */;
                long long time_elapsed_since_exit_request /* in seconds */, time_elapsed_since_pause_request /* in seconds */;
                std::chrono::milliseconds flash_period(current_bgm ? 255 : 170 /* a BGM specific value figured out by rough testing */);
                std::chrono::steady_clock::time_point game_begin_time = std::chrono::steady_clock::now(), last_time = game_begin_time, exit_request_time, pause_request_time;
                bool bg_refresh = true, prepare_exit = false, prepare_pause = false /* they must not be true at the same time */, ignore_mouse_hover_on_pause = (page | RESUMED_FILTER) == RESUMED, refresh_if_not_really_exit = false;
                /** Awaiting for mouse action while puzzle not solved */
                while (true) {
                    std::chrono::steady_clock::time_point time_this_loop = std::chrono::steady_clock::now();
                    if (prepare_exit) {
                        time_elapsed_since_exit_request = std::chrono::duration_cast<std::chrono::seconds>(time_this_loop - exit_request_time).count();
                        if (time_elapsed_since_exit_request >= countdown_before_exit) {
                            page = EXIT_PAGE;
                            break;
                        }
                    }
                    else if (prepare_pause) {
                        time_elapsed_since_pause_request = std::chrono::duration_cast<std::chrono::seconds>(time_this_loop - pause_request_time).count();
                        if (time_elapsed_since_pause_request >= countdown_before_pause) {
                            page = SQUARE_MODE | HUE_MODE | PAUSED_PAGE;
                            break;
                        }
                    }
                    if (bg_refresh) {
                        if (refresh_if_not_really_exit) {
                            square::fill_game_region(this_game, &hue_play);
                            square::draw_outer_border(this_game, 3, false);
                        }
                        if (prepare_pause || prepare_exit) {
                            settextcolor(RED);
                            setbkmode(OPAQUE);
                        }
                        if (prepare_pause) {
                            long long empty_block_begin_x = this_game->prim_sqr.begin_x + (this_game->get_empty() % this_game->get_N()) * (this_game->prim_sqr.size / this_game->get_N()), empty_block_begin_y = this_game->prim_sqr.begin_y + (this_game->get_empty() / this_game->get_N()) * (this_game->prim_sqr.size / this_game->get_N());
                            long long text_height;
                            switch (this_game->get_N()) {
                            case 1:
                            case 2:
                            case 3:
                            case 4:
                                text_height = 100;
                                break;
                            case 5:
                                text_height = 50;
                                break;
                            default:
                                text_height = 30;
                            } //TODO More testing
                            setbkcolor(WHITE);
                            settextstyle(text_height, 0, L"Consolas", 0, 0, 400, false, false, false);
                            long long tmp_x = empty_block_begin_x + this_game->prim_sqr.size / this_game->get_N() / 2, tmp_y = empty_block_begin_y + this_game->prim_sqr.size / this_game->get_N() / 2;
                            outtextxy(tmp_x - custom::consolas_width(text_height) / 2, tmp_y - text_height / 2, std::to_wstring(countdown_before_pause - time_elapsed_since_pause_request /* a single digit */).c_str());
                            settextstyle(20, 0, L"Candara", 0, 0, 400, false, false, false);
                            outtextxy(tmp_x - 33, tmp_y + text_height / 2 + 3, L"- Paused -"); //TODO Documentate as temporary pausing
                            square::draw_outer_border(this_game, 3, false);
                        }
                        else if (prepare_exit) {
                            setlinecolor(WHITE);
                            setfillcolor(0xe2e2e2);
                            setbkcolor(0xe2e2e2);
                            fillrectangle(418, 425, 578, 625);
                            settextstyle(150, 0, L"Consolas", 0, 0, 400, false, false, false);
                            outtextxy(500 - custom::consolas_width(150) / 2, 425, std::to_wstring(countdown_before_exit - time_elapsed_since_exit_request /* a single digit */).c_str());
                            settextstyle(40, 0, L"Candara", 0, 0, 400, false, false, false);
                            outtextxy(500 - 75, 425 + 150 + 3, L"- Quitting -");
                        }
                        if (prepare_pause || prepare_exit) {
                            setbkmode(TRANSPARENT);
                        }
                        FlushBatchDraw();
                        bg_refresh = false;
                    }
                    else {
                        if (std::chrono::duration_cast<std::chrono::milliseconds>(time_this_loop - last_time) >= flash_period)
                            last_time = time_this_loop, bg_refresh = true;
                    }
                    if (!(prepare_pause || prepare_exit)) /* put outside of the above refresh cycle for precision of timer */ {
                        if (std::chrono::duration_cast<std::chrono::milliseconds>(time_this_loop - game_begin_time + this_game->get_timer()) < this_game->get_time_limit()) {
                            long long tmp;
                            show_timer(3, tmp = std::chrono::duration_cast<std::chrono::milliseconds>(time_this_loop - game_begin_time + this_game->get_timer()).count(), tmp * 10 > std::chrono::duration_cast<std::chrono::milliseconds>(this_game->get_time_limit()).count() * 9) /* start to show red text when exceeded 90% the allowed time */;
                        }
                        else {
                            show_timer(3, std::chrono::duration_cast<std::chrono::milliseconds>(this_game->get_time_limit()).count(), true);
                            sound_after_game(false);
                            this_game->finalize_timer();
                            page = SQUARE_MODE | HUE_MODE | RESULT_PAGE; /* failed */
                            break;
                        }
                    }
                    bool peeked = peekmessage(&msg, EX_MOUSE | EX_KEY);
                    if (peeked) {
                        if (msg.message == WM_KEYUP /* suppose the keyup after keydown is successfully detected */ && msg.vkcode == VK_ESCAPE) {
                            page = START_PAGE;
                            break;
                        }
                        if (msg.message == WM_MOUSEMOVE /* I believe this suffices */)
                            ignore_mouse_hover_on_pause = ignore_mouse_hover_on_pause && custom::is_in_circ(&msg, 110, 38, 21) /* it is ok to continue the flow using msg possibly again */;
                        if (prepare_pause) {
                            if (msg.message == WM_KEYUP && msg.vkcode == VK_SPACE || msg.message == WM_LBUTTONUP && custom::is_in_circ(&msg, 110, 38, 21)) {
                                page = SQUARE_MODE | HUE_MODE | PAUSED_PAGE;
                                break;
                            }
                            if (!custom::is_in_circ(&msg, 110, 38, 21)) {
                                /** Resume the game */
                                square::fill_block(this_game, &hue_play, this_game->get_empty() / this_game->get_N(), this_game->get_empty() % this_game->get_N(), true);
                                square::draw_outer_border(this_game, 3, false);
                                game_begin_time = time_this_loop, prepare_pause = false;
                                /* ---------------- */
                                continue;
                            }
                        }
                        else if (prepare_exit) {
                            if (msg.message == WM_LBUTTONUP && custom::is_in_rect(&msg, 23, 24, 27, 32)) {
                                page = START_PAGE;
                                break;
                            }
                            if (!custom::is_in_rect(&msg, 23, 24, 27, 32)) {
                                /** Resume the game */
                                square::fill_block(this_game, &hue_play, this_game->get_empty() / this_game->get_N(), this_game->get_empty() % this_game->get_N(), true);
                                square::draw_outer_border(this_game, 3, false);
                                refresh_if_not_really_exit = true, game_begin_time = time_this_loop, prepare_exit = false, bg_refresh = true /* refresh immediately within the next loop */;
                                /* ---------------- */
                                continue;
                            }
                        }
                        else {
                            if (!square::move_if_asked(this_game, msg, hue_play, 3)) {
                                prepare_exit = custom::is_in_rect(&msg, 23, 24, 27, 32);
                                prepare_pause = !prepare_exit && !ignore_mouse_hover_on_pause && (msg.message == WM_KEYDOWN && msg.vkcode == VK_SPACE || custom::is_in_circ(&msg, 110, 38, 21));
                                if (prepare_exit || prepare_pause) {
                                    this_game->increment_timer(std::chrono::duration_cast<std::chrono::milliseconds>(time_this_loop - game_begin_time));
                                }
                                if (prepare_exit)
                                    exit_request_time = time_this_loop, bg_refresh = true /* refresh immediately within the next loop */;
                                if (prepare_pause)
                                    pause_request_time = time_this_loop, bg_refresh = true /* refresh immediately within the next loop */;
                            }
                        }
                    }
                    if (this_game->check_solved()) {
                        sound_after_game(true);
                        page = SQUARE_MODE | HUE_MODE | RESULT_PAGE; /* winned */
                        break;
                    }
                }
                /** ------------------------------------------------- */
                continue;
            }
        }
    }
    //TODO EndBatchDraw();
    closegraph();

    //cv::Mat image = cv::imread("Files\\Images\\Puzzle_Qiuzhen\\0.jpg", cv::IMREAD_UNCHANGED);
    //cv::blur(image, image, cv::Size(20, 20));
    //if (!image.empty() /* successfully loaded */) {
    //    cv::imwrite("Files\\Images\\TEMP\\0.jpg", image);
    //    cv::namedWindow("Image", cv::WINDOW_NORMAL);
    //    imshow("Image", image);
    //}
    //cv::waitKey(0); //1 is also ok, if no this line then imgae is not shown

    /** The vertices of the primitive triangle (made close to an equilateral triangle) */
    /*prim_tri.vertices[0] = std::pair<long long, long long>(200, 480 - 360);
    prim_tri.vertices[1] = std::pair<long long, long long>(200, 480 + 360);
    prim_tri.vertices[2] = std::pair<long long, long long>(824, 480);*/
    /** ------------------------------------------------------------------------------ */
    //fill_entire_graph(graph_buffer, &hue_play);
    //game(TRIANGLE_MODE, 2, &hue_play);
    return 0;
}