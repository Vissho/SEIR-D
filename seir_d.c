#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const double alpha_E = 0.999;
const double alpha_I = 0.999;
const double kappa = 0.042;
const double rho = 0.952;
const double beta = 0.999;
const double mu = 0.0188;
const double c_isol = 0;
const double E0 = 99;
const double R0 = 24;

const double N0 = 2798170.0;
const double S0 = N0 - E0 - R0;
const double I0 = 0.0;
const double D0 = 0.0;

const double Gamma = 0.0;
const double tau = 2;
const int SIZE = 5;
int cnt = 0;

double c(double t)
{
    double a = t;
    return 1 + c_isol * (1 - 2 / 5 * a);
}

int func(double t, double Y[SIZE], double res[SIZE])
{
    double S = Y[0];
    double E = Y[1];
    double I = Y[2];
    double R = Y[3];
    double D = Y[4];
    double N = S + E + I + R + D;
    double c_t = c(t - tau);
    double dSdt
            = -c_t * (alpha_I * S * I / N + alpha_E * S * E / N) + Gamma * R;
    double dEdt = c_t * (alpha_I * S * I / N + alpha_E * S * E / N)
            - (kappa + rho) * E;
    double dIdt = kappa * E - beta * I - mu * I;
    double dRdt = beta * I + rho * E - Gamma * R;
    double dDdt = mu * I;
    res[0] = dSdt;
    res[1] = dEdt;
    res[2] = dIdt;
    res[3] = dRdt;
    res[4] = dDdt;
    return 0;
}

int add_data(double** result, double t, double Y[SIZE])
{
    cnt++;
    result[cnt - 1][0] = t;
    for (int i = 0; i < SIZE; i++)
        result[cnt - 1][i + 1] = Y[i];

    return 0;
}

double max_elem(double err[SIZE])
{
    double max = 0;
    for (int i = 0; i < SIZE; i++)
        if (max < err[i])
            max = err[i];

    return max;
}

int euler_method(
        double t_0,
        double Y_0[SIZE],
        double h,
        double t_max,
        double eps,
        double** result)
{
    double t = t_0;
    double Y[SIZE];
    for (int i = 0; i < SIZE; i++)
        Y[i] = Y_0[i];

    add_data(result, t, Y);

    while (t < t_max) {
        double Y_loc[SIZE];
        double Y_corr[SIZE];
        for (int i = 0; i < SIZE; i++) {
            Y_loc[i] = Y[i];
            Y_corr[i] = Y[i];
        }

        double res1[SIZE];
        func(t, Y, res1);
        for (int i = 0; i < SIZE; i++) {
            Y_loc[i] += h * res1[i];
        }

        double res2[SIZE];
        func(t + h, Y_loc, res2);
        for (int i = 0; i < SIZE; i++) {
            Y_corr[i] += (h / 2) * (res1[i] + res2[i]);
        }

        double err[SIZE];
        for (int i = 0; i < SIZE; i++) {
            err[i] = fabs(Y_corr[i] - Y_loc[i]);
        }

        double max_err = max_elem(err);
        if (max_err < eps) {
            t = t + h;
            for (int i = 0; i < SIZE; i++)
                Y[i] = Y_corr[i];

            add_data(result, t, Y);
        } else {
            h /= 2;
        }
    }

    return 0;
}

int main()
{
    double t_0 = 0.0;
    double t_max = 90.0;
    double h = 1;
    double eps = 0.01;
    double Y_0[] = {S0, E0, I0, R0, D0};

    int sizes = 50001;
    double** result = (double**)malloc(sizeof(double*) * sizes);
    for (int i = 0; i < sizes; i++)
        result[i] = (double*)malloc(sizeof(double) * 6);

    euler_method(t_0, Y_0, h, t_max, eps, result);

    double fp = 0, ip = 0;

    for (int i = 0; i < cnt; i++) {
        // printf("%lf\t%lf\n", t_max, h);
        fp = modf(result[i][0], &ip);
        if (fp == 0)
            printf("%.0lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",
                   result[i][0],
                   result[i][1],
                   result[i][2],
                   result[i][3],
                   result[i][4],
                   result[i][5]);
    }

    for (int i = 0; i < sizes; i++)
        free(result[i]);
    free(result);

    return 0;
}
