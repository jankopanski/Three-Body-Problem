#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <assert.h>

typedef struct {
    double x, y, z, vx, vy, vz, ax, ay, az, axp, ayp, azp;
} Particle;

const double E0 = 1.0;
const double eps_sqrt = 4.69041575982343e-08;
const int buffer_size = 256;

void update_coordinates(Particle *p, const int num, const double deltatime) {
    for (int i = 0; i < num; i++) {
        p[i].x += (p[i].vx + 0.5 * p[i].ax * deltatime) * deltatime;
        p[i].y += (p[i].vy + 0.5 * p[i].ay * deltatime) * deltatime;
        p[i].z += (p[i].vz + 0.5 * p[i].az * deltatime) * deltatime;
    }
}

double norm1(double x1, double x2, double y1, double y2, double z1, double z2) {
    return fmax((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2), 1e-20);
}

double compute_triple_potential(Particle pi, Particle pj, Particle pk) {
    double rij2 = norm1(pi.x, pj.x, pi.y, pj.y, pi.z, pj.z);
    double rik2 = norm1(pi.x, pk.x, pi.y, pk.y, pi.z, pk.z);
    double rjk2 = norm1(pj.x, pk.x, pj.y, pk.y, pj.z, pk.z);
    double rij = sqrt(rij2);
    double rik = sqrt(rik2);
    double rjk = sqrt(rjk2);
    double rijk = rij * rik * rjk;
//    printf("%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf\n", pi.x, pi.y, pi.z, pj.x, pj.y, pj.z, pk.x, pk.y, pk.z);
//    printf("%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf\n", rij2, rik2, rjk2, rij, rik, rjk, rijk);
    return E0 * (1.0 / pow(rijk, 3) + 0.375 * (-rij2 + rik2 + rjk2) * (rij2 - rik2 + rjk2) * (rij2 + rik2 - rjk2) / pow(rijk, 5));
}

Particle move_x(Particle p, double h) {
    p.x += h;
    return p;
}

Particle move_y(Particle p, double h) {
    p.y += h;
    return p;
}

Particle move_z(Particle p, double h) {
    p.z += h;
    return p;
}

double compute_h(double x) {
    const double min_val = 1e-10;
    return fabs(x) < min_val ? min_val : eps_sqrt * x;
}

void update_particle_acceleration(Particle *pi, Particle *pj, Particle *pk) {
    double hx = compute_h(pi->x);
    double hy = compute_h(pi->y);
    double hz = compute_h(pi->z);
    volatile double delta_x = (pi->x + hx) - (pi->x - hx);
    volatile double delta_y = (pi->y + hy) - (pi->y - hy);
    volatile double delta_z = (pi->z + hz) - (pi->z - hz);
//    printf("%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf\n", hx, hy, hz, delta_x, delta_y, delta_z);
//    printf("%.16lf %.16lf\n", compute_triple_potential(move_x(*pi, hx), *pj, *pk), compute_triple_potential(move_x(*pi, -hx), *pj, *pk));
    pi->ax -= 2.0 * (compute_triple_potential(move_x(*pi, hx), *pj, *pk) - compute_triple_potential(move_x(*pi, -hx), *pj, *pk)) / delta_x;
    pi->ay -= 2.0 * (compute_triple_potential(move_y(*pi, hy), *pj, *pk) - compute_triple_potential(move_y(*pi, -hy), *pj, *pk)) / delta_y;
    pi->az -= 2.0 * (compute_triple_potential(move_z(*pi, hz), *pj, *pk) - compute_triple_potential(move_z(*pi, -hz), *pj, *pk)) / delta_z;
}

void update_acceleration(Particle *p, const int num) {
    for (int i = 0; i < num; i++) {
        p[i].axp = p[i].ax;
        p[i].ayp = p[i].ay;
        p[i].azp = p[i].az;
        p[i].ax = 0;
        p[i].ay = 0;
        p[i].az = 0;
    }
    for (int i = 0; i < num; i++) {
        for (int j = i + 1; j < num; j++) {
            for (int k = j + 1; k < num; k++) {
                update_particle_acceleration(p + i, p + j, p + k);
                update_particle_acceleration(p + j, p + i, p + k);
                update_particle_acceleration(p + k, p + i, p + j);
            }
        }
    }
}

void update_velocity(Particle *p, const int num, const double deltatime) {
    for (int i = 0; i < num; i++) {
        p[i].vx += 0.5 * (p[i].axp + p[i].ax) * deltatime;
        p[i].vy += 0.5 * (p[i].ayp + p[i].ay) * deltatime;
        p[i].vz += 0.5 * (p[i].azp + p[i].az) * deltatime;
    }
}

void print_particles(Particle *p, int num, int step, char *filename) {
    char output_filename[256];
    char step_str[32];
    sprintf(step_str, "%d", step);
    strcpy(output_filename, filename);
    strcat(output_filename, "_");
    strcat(output_filename, step_str);
    strcat(output_filename, ".txt");
    FILE *output_stream = fopen(output_filename, "w");
    for (int i = 0; i < num; i++) {
        fprintf(output_stream, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf\n", p[i].x, p[i].y, p[i].z, p[i].vx, p[i].vy, p[i].vz);
    }
    fclose(output_stream);
}

void compute(Particle *p, const int num, const int stepcount, const double deltatime, bool verbose, char *filename) {
    update_acceleration(p, num);
//    for (int i = 0; i < 8; i++) {
//        printf("%d %.16lf %.16lf %.16lf\n", i, p[i].ax, p[i].ay, p[i].az);
//    }
    for (int s = 1; s <= stepcount; s++) {
        update_coordinates(p, num, deltatime);
        update_acceleration(p, num);
//        for (int i = 0; i < 4; i++) {
//            printf("%d %.16lf %.16lf %.16lf\n", i, p[i].ax, p[i].ay, p[i].az);
//        }
        update_velocity(p, num, deltatime);
        if (verbose) {
            print_particles(p, num, s, filename);
        }
    }
    if (!verbose) {
        print_particles(p, num, stepcount, filename);
    }
}

int main(int argc, char *argv[]) {
    int stepcount = 0;
    double deltatime = 0;
    bool verbose = false;
    int num_particles = 0;
    char *in_file = NULL, *out_file = NULL;
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
    }

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
    assert(num_particles);
    rewind(input_stream);
    Particle *particles = malloc(num_particles * sizeof(Particle));
    for (int i = 0; i < num_particles; i++) {
        fscanf(input_stream, "%lf %lf %lf %lf %lf %lf",
               &particles[i].x, &particles[i].y, &particles[i].z,
               &particles[i].vx, &particles[i].vy, &particles[i].vz);
    }
    compute(particles, num_particles, stepcount, deltatime, verbose, out_file);
    return 0;
}
