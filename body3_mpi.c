#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <assert.h>
#include <mpi.h>


typedef struct {
    double x, y, z;
} Point;
MPI_Datatype dt_point;

const double E0 = 1.0;
const double eps_sqrt = 4.69041575982343e-08;
const Point zero_point = {0};


inline int prev_rank(int rank, int num_processes) {
    return rank > 0 ? rank - 1 : num_processes - 1;
}

inline int next_rank(int rank, int num_processes) {
    return (rank + 1) % num_processes;
}

inline int get_num_particles_from_rank(int rank, int num_processes, int num_particles) {
    return num_particles / num_processes + (rank < num_particles % num_processes);
}

inline double norm1(double x1, double x2, double y1, double y2, double z1, double z2) {
    return fmax((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2), 1e-20);
}

double compute_triple_potential(Point ci, Point cj, Point ck) {
    double rij2 = norm1(ci.x, cj.x, ci.y, cj.y, ci.z, cj.z);
    double rik2 = norm1(ci.x, ck.x, ci.y, ck.y, ci.z, ck.z);
    double rjk2 = norm1(cj.x, ck.x, cj.y, ck.y, cj.z, ck.z);
    double rij = sqrt(rij2);
    double rik = sqrt(rik2);
    double rjk = sqrt(rjk2);
    double rijk = rij * rik * rjk;
    return E0 * (1.0 / pow(rijk, 3)
                 + 0.375 * (-rij2 + rik2 + rjk2) * (rij2 - rik2 + rjk2) * (rij2 + rik2 - rjk2) / pow(rijk, 5));
}

Point move_x(Point c, double h) {
    c.x += h;
    return c;
}

Point move_y(Point c, double h) {
    c.y += h;
    return c;
}

Point move_z(Point c, double h) {
    c.z += h;
    return c;
}

double compute_h(double x) {
    const double min_val = 1e-10;
    return fabs(x) < min_val ? min_val : eps_sqrt * x;
}

Point compute_triple_force(Point ci, Point cj, Point ck) {
    double hx = compute_h(ci.x);
    double hy = compute_h(ci.y);
    double hz = compute_h(ci.z);
    volatile double delta_x = (ci.x + hx) - (ci.x - hx);
    volatile double delta_y = (ci.y + hy) - (ci.y - hy);
    volatile double delta_z = (ci.z + hz) - (ci.z - hz);
    Point f;
    f.x = -2.0 * (compute_triple_potential(move_x(ci, hx), cj, ck)
                  - compute_triple_potential(move_x(ci, -hx), cj, ck)) / delta_x;
    f.y = -2.0 * (compute_triple_potential(move_y(ci, hy), cj, ck)
                  - compute_triple_potential(move_y(ci, -hy), cj, ck)) / delta_y;
    f.z = -2.0 * (compute_triple_potential(move_z(ci, hz), cj, ck)
                  - compute_triple_potential(move_z(ci, -hz), cj, ck)) / delta_z;
    return f;
}

void compute_three_same_sets_forces(const Point *c0, Point *f0, int b0, int num_processes, int num_particles) {
    int set_size = get_num_particles_from_rank(b0, num_processes, num_particles);
    for (int i = 0; i < set_size; i++) {
        for (int j = i + 1; j < set_size; j++) {
            for (int k = j + 1; k < set_size; k++) {
                Point force_ijk = compute_triple_force(c0[i], c0[j], c0[k]);
                Point force_jik = compute_triple_force(c0[j], c0[i], c0[k]);
                Point force_kij = compute_triple_force(c0[k], c0[i], c0[j]);
                f0[i].x += force_ijk.x;
                f0[i].y += force_ijk.y;
                f0[i].z += force_ijk.z;
                f0[j].x += force_jik.x;
                f0[j].y += force_jik.y;
                f0[j].z += force_jik.z;
                f0[k].x += force_kij.x;
                f0[k].y += force_kij.y;
                f0[k].z += force_kij.z;
            }
        }
    }
}

void compute_two_same_sets_forces(const Point *c0, const Point *c1, Point *f0, Point *f1,
                                  int b0, int b1, int num_processes, int num_particles) {
    int set_size[2] = {get_num_particles_from_rank(b0, num_processes, num_particles),
                       get_num_particles_from_rank(b1, num_processes, num_particles)};
    for (int i = 0; i < set_size[0]; i++) {
        for (int j = i + 1; j < set_size[0]; j++) {
            for (int k = 0; k < set_size[1]; k++) {
                Point force_ijk = compute_triple_force(c0[i], c0[j], c1[k]);
                Point force_jik = compute_triple_force(c0[j], c0[i], c1[k]);
                Point force_kij = compute_triple_force(c1[k], c0[i], c0[j]);
                f0[i].x += force_ijk.x;
                f0[i].y += force_ijk.y;
                f0[i].z += force_ijk.z;
                f0[j].x += force_jik.x;
                f0[j].y += force_jik.y;
                f0[j].z += force_jik.z;
                f1[k].x += force_kij.x;
                f1[k].y += force_kij.y;
                f1[k].z += force_kij.z;
            }
        }
    }
}

void compute_unique_sets_forces(const Point *c0, const Point *c1, const Point *c2, Point *f0, Point *f1, Point *f2,
                                int b0, int b1, int b2, int num_processes, int num_particles) {
    int set_size[3] = {get_num_particles_from_rank(b0, num_processes, num_particles),
                       get_num_particles_from_rank(b1, num_processes, num_particles),
                       get_num_particles_from_rank(b2, num_processes, num_particles)};
    for (int i = 0; i < set_size[0]; i++) {
        for (int j = 0; j < set_size[1]; j++) {
            for (int k = 0; k < set_size[2]; k++) {
                Point force_ijk = compute_triple_force(c0[i], c1[j], c2[k]);
                Point force_jik = compute_triple_force(c1[j], c0[i], c2[k]);
                Point force_kij = compute_triple_force(c2[k], c0[i], c1[j]);
                f0[i].x += force_ijk.x;
                f0[i].y += force_ijk.y;
                f0[i].z += force_ijk.z;
                f1[j].x += force_jik.x;
                f1[j].y += force_jik.y;
                f1[j].z += force_jik.z;
                f2[k].x += force_kij.x;
                f2[k].y += force_kij.y;
                f2[k].z += force_kij.z;
            }
        }
    }
}

void shift_right(Point *c, Point *f, int *b, int rank, int num_processes, int num_particles_per_process) {
    const int coordinates_tag = 1;
    const int forces_tag = 2;
    MPI_Sendrecv_replace(c, num_particles_per_process, dt_point, next_rank(rank, num_processes), coordinates_tag,
                         prev_rank(rank, num_processes), coordinates_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv_replace(f, num_particles_per_process, dt_point, next_rank(rank, num_processes), forces_tag,
                         prev_rank(rank, num_processes), forces_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    *b = prev_rank(*b, num_processes);
}

void compute_forces(Point **c, Point **f, int *b, int rank, int num_processes, int num_particles, int num_particles_per_process) {
    int i = 0;
    for (int steps = num_processes - 3; steps >= num_processes % 3; steps -= 3) {
        for (int move = 0; move < steps; move++) {
            if (move != 0 || steps != num_processes - 3) {
                shift_right(c[i], f[i], &b[i], rank, num_processes, num_particles_per_process);
            }
            else {
                compute_three_same_sets_forces(c[1], f[1], b[1], num_processes, num_particles);
                compute_two_same_sets_forces(c[1], c[2], f[1], f[2], b[1], b[2], num_processes, num_particles);
                compute_two_same_sets_forces(c[0], c[2], f[0], f[2], b[0], b[2], num_processes, num_particles);
            }
            if (steps == num_processes - 3) {
                compute_two_same_sets_forces(c[1], c[0], f[1], f[0], b[1], b[0], num_processes, num_particles);
            }
            compute_unique_sets_forces(c[0], c[1], c[2], f[0], f[1], f[2], b[0], b[1], b[2], num_processes, num_particles);
        }
        i = (i + 1) % 3;
    }
    if (num_processes % 3 == 0) {
        i = prev_rank(i, 3);
        shift_right(c[i], f[i], &b[i], rank, num_processes, num_particles_per_process);
        if (rank / (num_processes / 3) == 0) {
            compute_unique_sets_forces(c[0], c[1], c[2], f[0], f[1], f[2], b[0], b[1], b[2], num_processes, num_particles);
        }
    }
}

void allocate_forces(Point **f, int num_particles_per_process) {
    for (int i = 0; i < 3; i++) {
        f[i] = malloc(num_particles_per_process * sizeof(Point));
        memset(f[i], 0, num_particles_per_process * sizeof(Point));
    }
}

void copy_coordinates(Point **c, const Point *my_coordinates, const int *b, int num_particles_per_process) {
    for (int i = 0; i < 3; i++) {
        c[i] = malloc(num_particles_per_process * sizeof(Point));
    }
    memcpy(c[1], my_coordinates, num_particles_per_process * sizeof(Point));
    MPI_Request requests[4];
    MPI_Status statuses[4];
    MPI_Irecv(c[0], num_particles_per_process, dt_point, b[0], 0, MPI_COMM_WORLD, &requests[0]);
    MPI_Irecv(c[2], num_particles_per_process, dt_point, b[2], 0, MPI_COMM_WORLD, &requests[1]);
    MPI_Isend(c[1], num_particles_per_process, dt_point, b[0], 0, MPI_COMM_WORLD, &requests[2]);
    MPI_Isend(c[1], num_particles_per_process, dt_point, b[2], 0, MPI_COMM_WORLD, &requests[3]);
    MPI_Waitall(4, requests, statuses);
}

void free_points(Point **p) {
    for (int i = 0; i < 3; i++) {
        free(p[i]);
    }
}

void collect_forces(Point **f, const int *b, Point *a, int num_particles_per_process) {
    MPI_Request requests[6];
    MPI_Status statuses[6];
    Point *fh[3]; // forces helpers
    for (int i = 0; i < 3; i++) {
        fh[i] = malloc(num_particles_per_process * sizeof(Point));
        memset(fh[i], 0, num_particles_per_process * sizeof(Point));
        MPI_Irecv(fh[i], num_particles_per_process, dt_point, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &requests[i]);
        MPI_Isend(f[i], num_particles_per_process, dt_point, b[i], 0, MPI_COMM_WORLD, &requests[3 + i]);
    }
    MPI_Waitall(6, requests, statuses);
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < num_particles_per_process; i++) {
        a[i].x = fh[0][i].x + fh[1][i].x + fh[2][i].x;
        a[i].y = fh[0][i].y + fh[1][i].y + fh[2][i].y;
        a[i].z = fh[0][i].z + fh[1][i].z + fh[2][i].z;
    }
    free_points(fh);
}

void update_acceleration(const Point *coordinates, Point *acceleration, Point *prev_acceleration,
                         int rank, int num_processes, int num_particles, int my_num_particles,
                         int num_particles_per_process) {
    for (int i = 0; i < my_num_particles; i++) {
        prev_acceleration[i].x = acceleration[i].x;
        prev_acceleration[i].y = acceleration[i].y;
        prev_acceleration[i].z = acceleration[i].z;
        acceleration[i].x = 0;
        acceleration[i].y = 0;
        acceleration[i].z = 0;
    }
    int b[3] = {prev_rank(rank, num_processes), rank, next_rank(rank, num_processes)};
    Point *c[3], *f[3]; // copy of coordinates, forces all equal to 0
    allocate_forces(f, num_particles_per_process);
    copy_coordinates(c, coordinates, b, num_particles_per_process);
    compute_forces(c, f, b, rank, num_processes, num_particles, num_particles_per_process);
    collect_forces(f, b, acceleration, num_particles_per_process);
    free_points(c);
    free_points(f);
}

void update_coordinates(Point *c, const Point *v, const Point *a, double deltatime, int my_num_particles) {
    for (int i = 0; i < my_num_particles; i++) {
        c[i].x += (v[i].x + 0.5 * a[i].x * deltatime) * deltatime;
        c[i].y += (v[i].y + 0.5 * a[i].y * deltatime) * deltatime;
        c[i].z += (v[i].z + 0.5 * a[i].z * deltatime) * deltatime;
    }
}

void update_velocity(Point *v, Point *a, Point *aprev, double deltatime, int my_num_particles) {
    for (int i = 0; i < my_num_particles; i++) {
        v[i].x += 0.5 * (aprev[i].x + a[i].x) * deltatime;
        v[i].y += 0.5 * (aprev[i].y + a[i].y) * deltatime;
        v[i].z += 0.5 * (aprev[i].z + a[i].z) * deltatime;
    }
}

void gather_and_print(Point *coordinates, Point *velocity, char *filename, int step, int rank, int num_processes,
                      int num_particles, int num_particles_per_process, int even_num_particles) {
    Point *output_coordinates = NULL;
    Point *output_velocity = NULL;
    if (rank == 0) {
        output_coordinates = malloc(even_num_particles * sizeof(Point));
        output_velocity = malloc(even_num_particles * sizeof(Point));
    }

    MPI_Gather(coordinates, num_particles_per_process, dt_point,
               output_coordinates, num_particles_per_process, dt_point, 0, MPI_COMM_WORLD);
    MPI_Gather(velocity, num_particles_per_process, dt_point,
               output_velocity, num_particles_per_process, dt_point, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        char output_filename[256];
        char step_str[32];
        sprintf(step_str, "%d", step);
        strcpy(output_filename, filename);
        strcat(output_filename, "_");
        strcat(output_filename, step_str);
        strcat(output_filename, ".txt");
        FILE *output_stream = fopen(output_filename, "w");
        for (int i = 0; i < even_num_particles; i++) {
            if (num_particles == even_num_particles
                || i / num_particles_per_process < num_particles % num_processes
                || i % num_particles_per_process != num_particles_per_process - 1) {
                fprintf(output_stream, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf\n",
                        output_coordinates[i].x, output_coordinates[i].y, output_coordinates[i].z,
                        output_velocity[i].x, output_velocity[i].y, output_velocity[i].z);
            }
        }
        fclose(output_stream);
        free(output_coordinates);
        free(output_velocity);
    }
}

void simulation(Point *coordinates, Point *velocity, Point *acceleration, int stepcount, double deltatime,
                int rank, int num_processes, int num_particles, int my_num_particles, int num_particles_per_process,
                int even_num_particles, char *filename, bool verbose) {
    Point *prev_acceleration = malloc(num_particles_per_process * sizeof(Point));
    update_acceleration(coordinates, acceleration, prev_acceleration,
                        rank, num_processes, num_particles, my_num_particles, num_particles_per_process);
    for (int s = 1; s <= stepcount; s++) {
        update_coordinates(coordinates, velocity, acceleration, deltatime, my_num_particles);
        update_acceleration(coordinates, acceleration, prev_acceleration,
                            rank, num_processes, num_particles, my_num_particles, num_particles_per_process);
        update_velocity(velocity, acceleration, prev_acceleration, deltatime, my_num_particles);
        if (verbose) {
            gather_and_print(coordinates, velocity, filename, s, rank, num_processes,
            num_particles, num_particles_per_process, even_num_particles);
        }
    }
    if (!verbose) {
        gather_and_print(coordinates, velocity, filename, stepcount, rank, num_processes,
                         num_particles, num_particles_per_process, even_num_particles);
    }
    free(prev_acceleration);
}

int main(int argc, char *argv[]) {
    int stepcount = 0;
    double deltatime = 0;
    bool verbose = false;
    char *in_file = NULL, *out_file = NULL;
    int num_particles = 0, num_particles_per_process, my_num_particles, even_num_particles;
    int num_processes, rank;
    Point *coordinates, *velocity, *acceleration;
    Point *input_coordinates = NULL, *input_velocity = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Type_contiguous(3, MPI_DOUBLE, &dt_point);
    MPI_Type_commit(&dt_point);
    assert(num_processes >= 4);

    if (argc >= 5 && argc <= 6) {
        in_file = argv[1];
        out_file = argv[2];
        stepcount = atoi(argv[3]);
        deltatime = atof(argv[4]);
        if (argc == 6 && !strcmp(argv[5], "-v")) {
            verbose = true;
        }
    }
    else {
        fprintf(stderr, "Invalid number of arguments\n");
        exit(EXIT_FAILURE);
    }

    if (rank == 0) {
        char *line = NULL;
        size_t len = 0;
        ssize_t nread;
        FILE *input_stream = fopen(in_file, "r");
        if (input_stream == NULL) {
            perror("fopen");
            exit(EXIT_FAILURE);
        }
        while ((nread = getline(&line, &len, input_stream)) != -1) {
            num_particles++;
        }
        fclose(input_stream);
        free(line);
        assert(num_particles);
        assert(num_particles >= num_processes);
    }

    MPI_Bcast(&num_particles, 1, MPI_INT, 0, MPI_COMM_WORLD);
    num_particles_per_process = (num_particles + num_processes - 1) / num_processes;
    my_num_particles = get_num_particles_from_rank(rank, num_processes, num_particles);
    even_num_particles = num_particles_per_process * num_processes;
    coordinates = malloc(num_particles_per_process * sizeof(Point));
    velocity = malloc(num_particles_per_process * sizeof(Point));
    acceleration = malloc(num_particles_per_process * sizeof(Point));
    memset(acceleration, 0, num_particles_per_process * sizeof(Point));

    if (rank == 0) {
        FILE *input_stream = fopen(in_file, "r");
        if (input_stream == NULL) {
            perror("fopen");
            exit(EXIT_FAILURE);
        }
        input_coordinates = malloc(even_num_particles * sizeof(Point));
        input_velocity = malloc(even_num_particles * sizeof(Point));
        for (int i = 0; i < num_processes; i++) {
            int lines_to_read = num_particles_per_process
                                - (num_particles % num_processes && i >= num_particles % num_processes);
            for (int j = 0; j < lines_to_read; j++) {
                int index = i * num_particles_per_process + j;
                fscanf(input_stream, "%lf %lf %lf %lf %lf %lf",
                       &input_coordinates[index].x, &input_coordinates[index].y, &input_coordinates[index].z,
                       &input_velocity[index].x, &input_velocity[index].y, &input_velocity[index].z);
            }
            if (lines_to_read < num_particles_per_process) {
                int index = i * num_particles_per_process + lines_to_read;
                input_coordinates[index] = zero_point;
                input_velocity[index] = zero_point;
            }
        }
        fclose(input_stream);
    }

    MPI_Scatter(input_coordinates, num_particles_per_process, dt_point,
                coordinates, num_particles_per_process, dt_point, 0, MPI_COMM_WORLD);
    MPI_Scatter(input_velocity, num_particles_per_process, dt_point,
                velocity, num_particles_per_process, dt_point, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        free(input_coordinates);
        free(input_velocity);
    }

    simulation(coordinates, velocity, acceleration, stepcount, deltatime, rank, num_processes,
               num_particles, my_num_particles, num_particles_per_process, even_num_particles,
               out_file, verbose);

    free(coordinates);
    free(velocity);
    free(acceleration);
    MPI_Finalize();
    return 0;
}
