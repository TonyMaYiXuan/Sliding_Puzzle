#include <graphics.h>
#include <conio.h>
#include <cstdlib>
#include <string>
#include <random>
#include <algorithm>
#include <fstream>
//#include "opencv-4.6.0/modules/core/include/opencv2/core.hpp"
//#include "opencv-4.6.0/modules/imgproc/include/opencv2/imgproc.hpp"
#pragma comment(lib, "winmm.lib")

class Triangle {
public:
    std::pair<long long, long long> vertices[3] /* (x1, y1), (x2, y2), (x3, y3), (x4, y4) */;
    short type /* 0 for upward and 1 for downward (take the primitive triangle as upward) */, boundary /* a number between 0 and 7 (inclusive), treated as a binary number */;
    long long target /* a number between 0 and custom_pow(4, ?) - 1 (inclusive), to be changed when swapping */, direction /* a number between 0 and 2 (inclusive), to be changed when swapping */;
};

class Square {
public:
    long long begin_x, begin_y, size;
    long long target /* a number between 0 and N * N - 1 (inclusive), to be changed when swapping */;
};

/** Global variable declarations */
long long graph_size /* in pixels */;
long long shuffle_times;
Triangle prim_tri /* primitive (original non-splited) triangle */;
Square prim_sqr /* primitive (original non-splited) square */;
/** ---------------------------- */

void initialization() {
    graph_size = 960;
    prim_sqr.begin_x = 60;
    prim_sqr.begin_y = 60;
    prim_sqr.size = 840; /* will be modified by the program to be a multiple of N */
    /** The vertices of the primitive triangle (made close to an equilateral triangle) */
    prim_tri.vertices[0] = std::pair<long long, long long>(200, 480 - 360);
    prim_tri.vertices[1] = std::pair<long long, long long>(200, 480 + 360);
    prim_tri.vertices[2] = std::pair<long long, long long>(824, 480);
    /** ------------------------------------------------------------------------------ */
    prim_tri.type = 0;
    prim_tri.boundary = 0b111;
    prim_tri.direction = 0;
    shuffle_times = 70000; /* for TRIANGLE_MODE, use a number close to the square of the number of blocks */
}

long long custom_pow(long long a, long long b /* non-negative integer */) {
    long long ans = 1;
    while (b)
        ans *= a, --b;
    return ans;
}

long long custom_rand(long long M) /* return an integer between 0 and M - 1 (inclusive) */ {
    std::random_device rd1 /* a seed for the random number engine */;
    std::mt19937 rd2(rd1()) /* standard mersenne_twister_engine seeded with rd1() */;
    // Use either rd1 or rd2 as the seed to feed the distribution below
    std::uniform_int_distribution<long long> distrib(0, M - 1);
    return distrib(rd2);
}

long long color_hue(long long i, long long j) {
    return BGR(RGB(j * 256 / graph_size, i * 256 / graph_size, 0));
}

#define TRIANGLE_MODE 0
#define SQUARE_MODE 1
namespace sqrgame {
    void fill_block(DWORD* pMem, Square* square_list, long long N, long long block_i, long long block_j, bool empty) {
        for (long long i = 0; i != prim_sqr.size / N; ++i)
            for (long long j = 0; j != prim_sqr.size / N; ++j) {
                long long goal_block_i = (square_list + block_i * N + block_j)->target / N, goal_block_j = (square_list + block_i * N + block_j)->target % N;
                pMem[(prim_sqr.begin_y + block_i * (prim_sqr.size / N) + i) * graph_size + (prim_sqr.begin_x + block_j * (prim_sqr.size / N) + j)] = empty ? WHITE : color_hue(prim_sqr.begin_y + goal_block_i * (prim_sqr.size / N) + i, prim_sqr.begin_x + goal_block_j * (prim_sqr.size / N) + j);
            }
    }
    void draw_outer_border(DWORD* pMem, long long N, long long line_radius, bool nice_boundary = true) {
        for (long long i = prim_sqr.begin_y; i != prim_sqr.begin_y + prim_sqr.size; ++i) {
            for (long long j = -line_radius; j < line_radius; ++j) {
                pMem[i * graph_size + (prim_sqr.begin_x + j)] = BLACK; /* left edge */
                pMem[i * graph_size + (prim_sqr.begin_x + prim_sqr.size + j)] = BLACK; /* right edge */
            }
        }
        long long tmp = nice_boundary * line_radius;
        for (long long i = prim_sqr.begin_y - tmp; i != prim_sqr.begin_y + prim_sqr.size + tmp; ++i) {
            for (long long j = -line_radius; j < line_radius; ++j) {
                pMem[(prim_sqr.begin_x + j) * graph_size + i] = BLACK; /* top edge */
                pMem[(prim_sqr.begin_x + prim_sqr.size + j) * graph_size + i] = BLACK; /* bottom edge */
            }
        }
    }
    void fill_game_region(DWORD* pMem, Square* square_list, long long N, long long empty_i, long long empty_j) {
        for (long long i = 0; i != N; ++i)
            for (long long j = 0; j != N; ++j)
                fill_block(pMem, square_list, N, i, j, i == empty_i && j == empty_j);
    }
    void swap(Square* square1, Square* square2) {
        long long tmp = square1->target;
        square1->target = square2->target;
        square2->target = tmp;
    }
    bool check_win(Square* square_list, long long N) {
        bool winned = true;
        for (long long i = 0; winned && i != N * N; ++i)
            winned = winned && square_list[i].target == i;
        return winned;
    }
}

namespace trigame {
    void triangle_construction(DWORD* pMem, Triangle** triangle_list, long long N, long long layer, long long k) {
        triangle_list[layer][k].target = k;
        if (layer == N)
            return;
        std::pair<long long, long long> A = triangle_list[layer][k].vertices[0], B = triangle_list[layer][k].vertices[1], C = triangle_list[layer][k].vertices[2];
        std::pair<long long, long long> F(A.first + B.first >> 1, A.second + B.second >> 1), E(C.first + A.first >> 1, C.second + A.second >> 1), D(B.first + C.first >> 1, B.second + C.second >> 1);
        triangle_list[layer + 1][4 * k].vertices[0] = A;
        triangle_list[layer + 1][4 * k].vertices[1] = F;
        triangle_list[layer + 1][4 * k].vertices[2] = E;
        triangle_list[layer + 1][4 * k].direction = triangle_list[layer][k].direction;
        triangle_list[layer + 1][4 * k].type = triangle_list[layer][k].type;
        triangle_list[layer + 1][4 * k].boundary = triangle_list[layer][k].boundary & 0b110;
        triangle_construction(pMem, triangle_list, N, layer + 1, 4 * k);
        triangle_list[layer + 1][4 * k + 1].vertices[0] = F;
        triangle_list[layer + 1][4 * k + 1].vertices[1] = B;
        triangle_list[layer + 1][4 * k + 1].vertices[2] = D;
        triangle_list[layer + 1][4 * k + 1].direction = triangle_list[layer][k].direction;
        triangle_list[layer + 1][4 * k + 1].type = triangle_list[layer][k].type;
        triangle_list[layer + 1][4 * k + 1].boundary = triangle_list[layer][k].boundary & 0b101;
        triangle_construction(pMem, triangle_list, N, layer + 1, 4 * k + 1);
        triangle_list[layer + 1][4 * k + 2].vertices[0] = E;
        triangle_list[layer + 1][4 * k + 2].vertices[1] = D;
        triangle_list[layer + 1][4 * k + 2].vertices[2] = C;
        triangle_list[layer + 1][4 * k + 2].direction = triangle_list[layer][k].direction;
        triangle_list[layer + 1][4 * k + 2].type = triangle_list[layer][k].type;
        triangle_list[layer + 1][4 * k + 2].boundary = triangle_list[layer][k].boundary & 0b011;
        triangle_construction(pMem, triangle_list, N, layer + 1, 4 * k + 2);
        triangle_list[layer + 1][4 * k + 3].vertices[0] = D;
        triangle_list[layer + 1][4 * k + 3].vertices[1] = E;
        triangle_list[layer + 1][4 * k + 3].vertices[2] = F;
        triangle_list[layer + 1][4 * k + 3].direction = triangle_list[layer][k].direction;
        triangle_list[layer + 1][4 * k + 3].type = (triangle_list[layer][k].type + 1) % 2;
        triangle_list[layer + 1][4 * k + 3].boundary = 0;
        triangle_construction(pMem, triangle_list, N, layer + 1, 4 * k + 3);
    }
    bool comp_x(const std::pair<long long, long long> A, const std::pair<long long, long long> B) {
        return (A.first < B.first) || (A.first == B.first && A.second < B.second);
    }
    bool comp_y(const std::pair<long long, long long> A, const std::pair<long long, long long> B) {
        return (A.second < B.second) || (A.second == B.second && A.first < B.first);
    }
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
    void fill_block(DWORD* pMem, Triangle** triangle_list, long long* pixel_status, long long N, long long layer, long long k, short mode) {
        /** Make a copy of vertices and sort them in custom order to find min / max */
        std::pair<long long, long long>* vertices = new std::pair<long long, long long>[3];
        for (long long i = 0; i != 3; ++i)
            vertices[i] = triangle_list[layer][k].vertices[i];
        std::sort(vertices, vertices + 3, comp_x);
        long long j_begin = vertices[0].first, j_end = vertices[2].first;
        std::sort(vertices, vertices + 3, comp_y);
        /** ----------------------------------------------------------------------- */
        for (long long i = vertices[0].second; i != vertices[2].second; ++i) {
            for (long long j = j_begin; j != j_end; ++j) {
                long long cnt = 0;
                cnt += ray_in_segment(std::pair<long long, long long>(j, i), vertices[0], vertices[1]);
                cnt += ray_in_segment(std::pair<long long, long long>(j, i), vertices[1], vertices[2]);
                cnt += ray_in_segment(std::pair<long long, long long>(j, i), vertices[2], vertices[0]);
                if (cnt == 1) {
                    /** The center of pixel (i, j) lies in the interior of the triangle */
                    if (mode & ASSIGN_PIXEL_STATUS)
                        pixel_status[i * graph_size + j] = k;
                    else if (mode & LEAVE_EMPTY)
                        pMem[i * graph_size + j] = WHITE;
                    else {
                        long long target_k = triangle_list[layer][k].target;
                        long long tmp_dir = triangle_list[layer][k].direction;
                        if (triangle_list[layer][k].type != triangle_list[layer][target_k].type)
                            triangle_list[layer][k].vertices[tmp_dir] = std::pair<long long, long long>(triangle_list[layer][k].vertices[(tmp_dir + 1) % 3].first + triangle_list[layer][k].vertices[(tmp_dir + 2) % 3].first - triangle_list[layer][k].vertices[tmp_dir].first, triangle_list[layer][k].vertices[(tmp_dir + 1) % 3].second + triangle_list[layer][k].vertices[(tmp_dir + 2) % 3].second - triangle_list[layer][k].vertices[tmp_dir].second);
                        std::pair<long double, long double> k_centroid((triangle_list[layer][k].vertices[0].first + triangle_list[layer][k].vertices[1].first + triangle_list[layer][k].vertices[2].first) / 3.0L, (triangle_list[layer][k].vertices[0].second + triangle_list[layer][k].vertices[1].second + triangle_list[layer][k].vertices[2].second) / 3.0F); /* reflected if necessary */
                        std::pair<long double, long double> target_k_centroid((triangle_list[layer][target_k].vertices[0].first + triangle_list[layer][target_k].vertices[1].first + triangle_list[layer][target_k].vertices[2].first) / 3.0L, (triangle_list[layer][target_k].vertices[0].second + triangle_list[layer][target_k].vertices[1].second + triangle_list[layer][target_k].vertices[2].second) / 3.0F);
                        long double rotate_angle = 3.1415926535897932384626433832795L * 4 * tmp_dir / 3;
                        std::pair<long double, long double> point(j + 0.5L, i + 0.5L);
                        if (triangle_list[layer][k].type != triangle_list[layer][target_k].type)
                            point = std::pair<long double, long double>(triangle_list[layer][k].vertices[(tmp_dir + 1) % 3].first + triangle_list[layer][k].vertices[(tmp_dir + 2) % 3].first - point.first, triangle_list[layer][k].vertices[(tmp_dir + 1) % 3].second + triangle_list[layer][k].vertices[(tmp_dir + 2) % 3].second - point.second);
                        std::pair<long long, long long> target_point = vector_minus(vector_add(rotate_point(vector_minus(point, k_centroid), rotate_angle), target_k_centroid), std::pair<long double, long double>(0.5L, 0.5L));
                        if (triangle_list[layer][k].type != triangle_list[layer][target_k].type)
                            triangle_list[layer][k].vertices[tmp_dir] = std::pair<long long, long long>(triangle_list[layer][k].vertices[(tmp_dir + 1) % 3].first + triangle_list[layer][k].vertices[(tmp_dir + 2) % 3].first - triangle_list[layer][k].vertices[tmp_dir].first, triangle_list[layer][k].vertices[(tmp_dir + 1) % 3].second + triangle_list[layer][k].vertices[(tmp_dir + 2) % 3].second - triangle_list[layer][k].vertices[tmp_dir].second);
                        pMem[i * graph_size + j] = color_hue(target_point.second, target_point.first);
                    }
                    /** --------------------------------------------------------------- */
                }
            }
        }
    }
    const long long ALGO_DIST_TO_SIDE = 0, ALGO_SMALLER_AND_BIGGER_INTERIOR = 1; /* similar to a macro */
    void draw_outer_border(DWORD* pMem, long long N, long long line_radius /* distance is calculated to the center of the pixel instead of the point in the pixel farthest away from the line */, short algo = ALGO_SMALLER_AND_BIGGER_INTERIOR) {
        std::pair<long long, long long> vertices[3];
        for (long long i = 0; i != 3; ++i)
            vertices[i] = prim_tri.vertices[i];
        if (algo == ALGO_SMALLER_AND_BIGGER_INTERIOR) {
            std::pair<long double, long double> vertices_pm[2][3];
            long long j_begin[2], i_begin[2], j_end[2], i_end[2];
            for (short k = 0; k != 3; ++k) {
                long double x = vertices[k].first - (vertices[0].first + vertices[1].first + vertices[2].first) / 3.0L, y = vertices[k].second - (vertices[0].second + vertices[1].second + vertices[2].second) / 3.0L;
                std::pair<long double, long double> delta(line_radius * x * sqrtl(3.0L / (x * x + y * y)), line_radius * y * sqrtl(3.0L / (x * x + y * y)));
                vertices_pm[0][k] = vector_add(vertices[k], delta), vertices_pm[1][k] = vector_minus(vertices[k], delta);
            }
            for (short k = 0; k != 2; ++k)
                j_begin[k] = min_element(vertices_pm[k], vertices_pm[k] + 3, comp_x)->first, j_end[k] = max_element(vertices_pm[k], vertices_pm[k] + 3, comp_x)->first, i_begin[k] = min_element(vertices_pm[k], vertices_pm[k] + 3, comp_y)->second, i_end[k] = max_element(vertices_pm[k], vertices_pm[k] + 3, comp_y)->second;
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
                        pMem[i * graph_size + j] = BLACK;
                }
        }
        if (algo == ALGO_DIST_TO_SIDE) {
            for (long long i = 0; i != graph_size; ++i)
                for (long long j = 0; j != graph_size; ++j) {
                    long double dist_sidelines[3];
                    short K = 0, k = 0;
                    do {
                        bool flag;
                        /** Calculate whether the center of the pixel is near the side on top of being close to the extended sideline (details omitted) */
                        long long Ax = 2 * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * vertices[(k + 1) % 3].first + (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (2 * j + 1) + (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (2 * vertices[(k + 1) % 3].second - 2 * i - 1), Ay = (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (2 * i + 1) + 2 * (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * vertices[(k + 1) % 3].second + (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (2 * vertices[(k + 1) % 3].first - 2 * j - 1), Bx = 2 * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * vertices[(k + 2) % 3].first + (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (2 * j + 1) + (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (2 * vertices[(k + 2) % 3].second - 2 * i - 1), By = (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (2 * i + 1) + 2 * (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * vertices[(k + 2) % 3].second + (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (2 * vertices[(k + 2) % 3].first - 2 * j - 1), delta = (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) + (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second);
                        if (Ax == Bx)
                            flag = (Ay <= delta * (2 * i + 1) && delta * (2 * i + 1) <= By) || (Ay >= delta * (2 * i + 1) && delta * (2 * i + 1) >= By);
                        else
                            flag = (Ax <= delta * (2 * j + 1) && delta * (2 * j + 1) <= Bx) || (Ax >= delta * (2 * j + 1) && delta * (2 * j + 1) >= Bx);
                        /** --------------------------------------------------------------------------------------------------------------------------- */
                        if (flag)
                            dist_sidelines[K++] = abs((vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (j + 0.5L) - (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (i + 0.5L) + (vertices[(k + 2) % 3].first * vertices[(k + 1) % 3].second - vertices[(k + 1) % 3].first * vertices[(k + 2) % 3].second)) / sqrtl((vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) * (vertices[(k + 2) % 3].first - vertices[(k + 1) % 3].first) + (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second) * (vertices[(k + 2) % 3].second - vertices[(k + 1) % 3].second)) /* make sure it is floating point division */;
                    } while (++k != 3);
                    if (K && *std::min_element(dist_sidelines, dist_sidelines + K) < line_radius)
                        pMem[i * graph_size + j] = BLACK;
                }
        }

    }

    void fill_game_region(DWORD* pMem, Triangle** triangle_list, long long *pixel_status, long long N, bool assign_pixel_status, long long empty = -1 /* used only if assign_pixel_status is false */) {
        for (long long i = 0; i != custom_pow(4, N - 1); ++i)
            fill_block(pMem, triangle_list, pixel_status, N, N - 1, i, assign_pixel_status? ASSIGN_PIXEL_STATUS: (i == empty? LEAVE_EMPTY: 0));
    }
    void swap(Triangle* triangle1, Triangle* triangle2, long long adjacent_type) {
        Triangle* tmp_triangle = new Triangle;
        *tmp_triangle = *triangle1;
        triangle1->target = triangle2->target;
        triangle1->direction = (6 - adjacent_type - triangle2->direction) % 3;
        triangle2->target = tmp_triangle->target;
        triangle2->direction = (6 - adjacent_type - tmp_triangle->direction) % 3;
    }
    bool check_adjacent_and_swap(Triangle** triangle_list, long long N, long long target, long long empty) {
        bool adjacent = true;
        if (triangle_list[N - 1][empty].vertices[1] == triangle_list[N - 1][target].vertices[2] && triangle_list[N - 1][empty].vertices[2] == triangle_list[N - 1][target].vertices[1])
            trigame::swap(*(triangle_list + N - 1) + empty, *(triangle_list + N - 1) + target, 0);
        else if (triangle_list[N - 1][empty].vertices[2] == triangle_list[N - 1][target].vertices[0] && triangle_list[N - 1][empty].vertices[0] == triangle_list[N - 1][target].vertices[2])
            trigame::swap(*(triangle_list + N - 1) + empty, *(triangle_list + N - 1) + target, 1);
        else if (triangle_list[N - 1][empty].vertices[0] == triangle_list[N - 1][target].vertices[1] && triangle_list[N - 1][empty].vertices[1] == triangle_list[N - 1][target].vertices[0])
            trigame::swap(*(triangle_list + N - 1) + empty, *(triangle_list + N - 1) + target, 2);
        else
            adjacent = false;
        return adjacent;
    }
    bool check_win(Triangle** triangle_list, long long N) {
        bool winned = true;
        for (long long i = 0; winned && i != custom_pow(4, N - 1); ++i)
            winned = winned && triangle_list[N - 1][i].target == i;
        return winned;
    }
}

/** ----- ARCHIVED -----
void write_err_log(std::string err) {
    std::ofstream err_log;
    err_log.open("Files\\ERR_LOG.txt", std::fstream::app);
    err_log << err;
    err_log.close();
}
*/

void game(short mode, long long N /* positive integer, meaningless practically if equals 1 */) {
    prim_sqr.size = prim_sqr.size / N * N; //TODO maybe changed quite strangely when such permanant changes was applied for many times
    initgraph(graph_size, graph_size);
    DWORD* pMem = GetImageBuffer();
    /** Fill the entire graph by a hue of color */
    for (long long i = 0; i != graph_size; ++i) {
        for (long long j = 0; j != graph_size; ++j) {
            pMem[i * graph_size + j] = color_hue(i, j); //TODO other color hue mode, but compatible with white texts
        }
    }
    /** --------------------------------------- */
    if (mode == SQUARE_MODE) {
        Square* square_list = new Square[N * N];
        long long block_size = prim_sqr.size / N;
        for (long long i = 0; i != N; ++i) {
            for (long long j = 0; j != N; ++j) {
                square_list[i * N + j].begin_x = prim_sqr.begin_x + j * block_size;
                square_list[i * N + j].begin_y = prim_sqr.begin_y + i * block_size;
                square_list[i * N + j].size = block_size /* maybe not used */;
                square_list[i * N + j].target = i * N + j;
            }
        }
        long long tmp_empty = custom_rand(N * N);
        long long empty_i = tmp_empty / N, empty_j = tmp_empty % N;
        /** Shuffle the game region */
        for (long long k = 0; k != shuffle_times; ++k) { //TODO Need to consider the case where no shuffle is done
            long long tmp_rand = custom_rand(4), delta_i, delta_j;
            switch (tmp_rand) {
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
            if (delta_i + empty_i >= 0 && delta_i + empty_i != N && delta_j + empty_j >= 0 && delta_j + empty_j != N) {
                sqrgame::swap(square_list + empty_i * N + empty_j, square_list + (empty_i + delta_i) * N + (empty_j + delta_j));
                empty_i += delta_i, empty_j += delta_j;
            }
        }
        /** ----------------------- */
        sqrgame::fill_game_region(pMem, square_list, N, empty_i, empty_j);
        sqrgame::draw_outer_border(pMem, N, 3);
        //FlushBatchDraw();//TODO Learn where it should be put at
        /** Awaiting for mouse action while not winned */
        ExMessage msg;
        while (!sqrgame::check_win(square_list, N)) {
            peekmessage(&msg, EX_MOUSE); //TODO peekmessage returns a boolean value, need to use or not?
            if (msg.message == WM_LBUTTONDOWN) {
                if (msg.x >= prim_sqr.begin_x && msg.x < prim_sqr.begin_x + prim_sqr.size && msg.y >= prim_sqr.begin_y && msg.y < prim_sqr.begin_y + prim_sqr.size) {
                    long long new_i = (msg.y - prim_sqr.begin_y) / block_size, new_j = (msg.x - prim_sqr.begin_x) / block_size, delta_i = new_i - empty_i, delta_j = new_j - empty_j;
                    if (delta_i * delta_i + delta_j * delta_j == 1) {
                        sqrgame::swap(square_list + empty_i * N + empty_j, square_list + new_i * N + new_j);
                        sqrgame::fill_block(pMem, square_list, N, empty_i, empty_j, false);
                        empty_i += delta_i, empty_j += delta_j;
                        //TODO Sound effects
                        sqrgame::fill_block(pMem, square_list, N, empty_i, empty_j, true);
                        sqrgame::draw_outer_border(pMem, N, 3, false /* not necessary */);
                    }
                }
            }
        }
        /** ------------------------------------------ */
        //TODO Wait for some time after winned
        sqrgame::fill_block(pMem, square_list, N, empty_i, empty_j, false);
        sqrgame::draw_outer_border(pMem, N, 3, false /* not necessary */);
        //TODO Sound effects for winned
    }
    else if (mode == TRIANGLE_MODE) {
        Triangle** triangle_list = new Triangle * [N];
        long long* pixel_status = new long long [graph_size * graph_size];
        for (long long i = 0; i != N; ++i)
            triangle_list[i] = new Triangle[custom_pow(4, i)];
        for (long long i = 0; i != graph_size * graph_size; ++i)
            pixel_status[i] = -1;
        **triangle_list = prim_tri;
        trigame::triangle_construction(pMem, triangle_list, N - 1, 0, 0);
        trigame::fill_game_region(pMem, triangle_list, pixel_status, N, true);
        long long empty = custom_rand(custom_pow(4, N - 1));
        /** Shuffle the game region */
        for (long long k = 0; k != shuffle_times; ++k) { //TODO Need to consider the case where no shuffle is done
            long long tmp_rand = custom_rand(3);
            std::pair<long double, long double> empty_centroid((triangle_list[N - 1][empty].vertices[0].first + triangle_list[N - 1][empty].vertices[1].first + triangle_list[N - 1][empty].vertices[2].first) / 3.0L, (triangle_list[N - 1][empty].vertices[0].second + triangle_list[N - 1][empty].vertices[1].second + triangle_list[N - 1][empty].vertices[2].second) / 3.0F);
            std::pair<long long, long long> target_centroid(2.0L * empty_centroid.first - triangle_list[N - 1][empty].vertices[tmp_rand].first - 0.5L, 2.0L * empty_centroid.second - triangle_list[N - 1][empty].vertices[tmp_rand].second - 0.5L);
            //TODO the well behavior of this routine depends somehow on the precision
            long long target;
            if (target_centroid.first >= 0 && target_centroid.first < graph_size && target_centroid.second >= 0 && target_centroid.second < graph_size && (target = pixel_status[target_centroid.second * graph_size + target_centroid.first]) != -1) {
                if (target == empty /* I believe this should not happen */)
                    throw std::runtime_error("Should not equal. empty = " + std::to_string(empty) + "; target = " + std::to_string(target) + ".\n");
                else if (!trigame::check_adjacent_and_swap(triangle_list, N, target, empty)) /* I believe the line below should not execute */ {
                    throw std::runtime_error("Not adjacent. empty = " + std::to_string(empty) + "; target = " + std::to_string(target));
                }
                empty = target;
            }
        }
        /** ----------------------- */
        trigame::fill_game_region(pMem, triangle_list, pixel_status, N, false, empty);
        trigame::draw_outer_border(pMem, N, 2);
        /** Awaiting for mouse action while not winned */
        ExMessage msg;
        while (!trigame::check_win(triangle_list, N)) {
            peekmessage(&msg, EX_MOUSE);
            if (msg.message == WM_LBUTTONDOWN) {
                long long target;
                if ((target = pixel_status[msg.y * graph_size + msg.x]) != -1 && target != empty && trigame::check_adjacent_and_swap(triangle_list, N, target, empty)) {
                    trigame::fill_block(pMem, triangle_list, pixel_status, N, N - 1, empty, 0);
                    bool need_draw_border = triangle_list[N - 1][empty].boundary > 0;
                    empty = target;
                    //TODO Sound effects
                    trigame::fill_block(pMem, triangle_list, pixel_status, N, N - 1, empty, trigame::LEAVE_EMPTY);
                    if (need_draw_border)
                        trigame::draw_outer_border(pMem, N, 2);
                }
            }
        }
        /** ------------------------------------------ */
        //TODO Wait for some time after winned
        trigame::fill_block(pMem, triangle_list, pixel_status, N, N - 1, empty, 0);
        if (triangle_list[N - 1][empty].boundary > 0)
            trigame::draw_outer_border(pMem, N, 2);
        //TODO Sound effects for winned
    }
}

int main() {
    initialization();
    WCHAR path[MAX_PATH];
    GetCurrentDirectoryW(MAX_PATH, path);
    std::wstring wstrpath(path);
    wstrpath += L"\\Files\\BGM_Nova_Ahrix.wav";
    PlaySoundW(wstrpath.c_str(), NULL, SND_ASYNC | SND_LOOP); //only runs after game ends //TODO need to avoid absolute path but this function needs absolute path? see Lv Xinyao's project
    game(TRIANGLE_MODE, 4);
    // TRIANGLE_MODE 1 means 1 layer
    _getch(); //TODO Put elsewhere
    closegraph(); //TODO Put elsewhere
    return 0;
}
