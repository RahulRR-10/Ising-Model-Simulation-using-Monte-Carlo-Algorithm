#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define kB 1.0
#define NUM_TEMPERATURES 100   
#define TEMP_START 0.0          
#define TEMP_END 300.0         
#define TEMP_STEP ((TEMP_END - TEMP_START) / NUM_TEMPERATURES)
#define L 200                   
#define J 1.5                   
#define MCS 50000               

void monteCarloStep(int *spin, double temperature);
double calculateAverageMagnetization(int *spin);

int main() {
    srand(time(NULL));

    FILE *dataFile = fopen("magnetization_vs_temperature.dat", "w");
    if (!dataFile) {
        printf("Error opening data file.\n");
        return 1;
    }

    for (double temperature = TEMP_START; temperature <= TEMP_END; temperature += TEMP_STEP) {
        int *spin = (int *)malloc(L * sizeof(int));
        if (spin == NULL) {
            printf("Memory allocation failed.\n");
            return 1;
        }

        for (int i = 0; i < L; i++) {
            spin[i] = (rand() % 2) * 2 - 1; 
        }

        for (int step = 0; step < MCS; step++) {
            monteCarloStep(spin, temperature);
        }

        double avg_magnetization = calculateAverageMagnetization(spin);

        fprintf(dataFile, "%.4f %.4f\n", temperature, avg_magnetization);

        free(spin);
    }

    fclose(dataFile);

    FILE *gnuplotPipe = popen("gnuplot -persist", "w");
    if (!gnuplotPipe) {
        printf("Error opening Gnuplot.\n");
        return 1;
    }

    fprintf(gnuplotPipe, "set xlabel 'Temperature (K)'\n");
    fprintf(gnuplotPipe, "set ylabel 'Average Magnetization'\n");
    fprintf(gnuplotPipe, "plot 'magnetization_vs_temperature.dat' with lines title 'Average Magnetization vs Temperature'\n");

    printf("Press Enter to exit...");
    getchar();

    fprintf(gnuplotPipe, "exit\n");
    pclose(gnuplotPipe);

    return 0;
}

void monteCarloStep(int *spin, double temperature) {
    for (int i = 0; i < L; i++) {
        int random_site = rand() % L;
        int delta_E = 2 * J * spin[random_site] * (spin[(random_site + 1) % L] + spin[(random_site - 1 + L) % L]);
        double acceptance_probability = exp(-delta_E / (kB * temperature));
        
        if (delta_E < 0 || (double)rand() / RAND_MAX < acceptance_probability) {
            spin[random_site] *= -1;
        }
    }
}

double calculateAverageMagnetization(int *spin) {
    double mag = 0.0;
    for (int i = 0; i < L; i++) {
        mag += (double)spin[i] / L;
    }
    return mag;
}
