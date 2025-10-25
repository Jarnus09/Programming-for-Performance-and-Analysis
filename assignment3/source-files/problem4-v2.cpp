#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define NSEC_SEC_MUL (1.0e9)

void gridloopsearch(
    double dd1, double dd2, double dd3, double dd4, double dd5, double dd6, double dd7, double dd8,
    double dd9, double dd10, double dd11, double dd12, double dd13, double dd14, double dd15,
    double dd16, double dd17, double dd18, double dd19, double dd20, double dd21, double dd22,
    double dd23, double dd24, double dd25, double dd26, double dd27, double dd28, double dd29,
    double dd30, double c11, double c12, double c13, double c14, double c15, double c16, double c17,
    double c18, double c19, double c110, double d1, double ey1, double c21, double c22, double c23,
    double c24, double c25, double c26, double c27, double c28, double c29, double c210, double d2,
    double ey2, double c31, double c32, double c33, double c34, double c35, double c36, double c37,
    double c38, double c39, double c310, double d3, double ey3, double c41, double c42, double c43,
    double c44, double c45, double c46, double c47, double c48, double c49, double c410, double d4,
    double ey4, double c51, double c52, double c53, double c54, double c55, double c56, double c57,
    double c58, double c59, double c510, double d5, double ey5, double c61, double c62, double c63,
    double c64, double c65, double c66, double c67, double c68, double c69, double c610, double d6,
    double ey6, double c71, double c72, double c73, double c74, double c75, double c76, double c77,
    double c78, double c79, double c710, double d7, double ey7, double c81, double c82, double c83,
    double c84, double c85, double c86, double c87, double c88, double c89, double c810, double d8,
    double ey8, double c91, double c92, double c93, double c94, double c95, double c96, double c97,
    double c98, double c99, double c910, double d9, double ey9, double c101, double c102,
    double c103, double c104, double c105, double c106, double c107, double c108, double c109,
    double c1010, double d10, double ey10, double kk);

struct timespec begin_grid, end_main;

// to store values of disp.txt
double a[120];

// to store values of grid.txt
double b[30];

int main()
{
    int i, j;

    i = 0;
    FILE *fp = fopen("./disp.txt", "r");
    if (fp == NULL)
    {
        printf("Error: could not open file\n");
        return 1;
    }

    while (!feof(fp))
    {
        if (!fscanf(fp, "%lf", &a[i]))
        {
            printf("Error: fscanf failed while reading disp.txt\n");
            exit(EXIT_FAILURE);
        }
        i++;
    }
    fclose(fp);

    // read grid file
    j = 0;
    FILE *fpq = fopen("./grid.txt", "r");
    if (fpq == NULL)
    {
        printf("Error: could not open file\n");
        return 1;
    }

    while (!feof(fpq))
    {
        if (!fscanf(fpq, "%lf", &b[j]))
        {
            printf("Error: fscanf failed while reading grid.txt\n");
            exit(EXIT_FAILURE);
        }
        j++;
    }
    fclose(fpq);

    // grid value initialize
    // initialize value of kk;
    double kk = 0.3;

    clock_gettime(CLOCK_MONOTONIC_RAW, &begin_grid);
    gridloopsearch(b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9], b[10], b[11], b[12],
                   b[13], b[14], b[15], b[16], b[17], b[18], b[19], b[20], b[21], b[22], b[23], b[24],
                   b[25], b[26], b[27], b[28], b[29], a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7],
                   a[8], a[9], a[10], a[11], a[12], a[13], a[14], a[15], a[16], a[17], a[18], a[19],
                   a[20], a[21], a[22], a[23], a[24], a[25], a[26], a[27], a[28], a[29], a[30], a[31],
                   a[32], a[33], a[34], a[35], a[36], a[37], a[38], a[39], a[40], a[41], a[42], a[43],
                   a[44], a[45], a[46], a[47], a[48], a[49], a[50], a[51], a[52], a[53], a[54], a[55],
                   a[56], a[57], a[58], a[59], a[60], a[61], a[62], a[63], a[64], a[65], a[66], a[67],
                   a[68], a[69], a[70], a[71], a[72], a[73], a[74], a[75], a[76], a[77], a[78], a[79],
                   a[80], a[81], a[82], a[83], a[84], a[85], a[86], a[87], a[88], a[89], a[90], a[91],
                   a[92], a[93], a[94], a[95], a[96], a[97], a[98], a[99], a[100], a[101], a[102],
                   a[103], a[104], a[105], a[106], a[107], a[108], a[109], a[110], a[111], a[112],
                   a[113], a[114], a[115], a[116], a[117], a[118], a[119], kk);
    clock_gettime(CLOCK_MONOTONIC_RAW, &end_main);
    printf("Total time = %f seconds\n", (end_main.tv_nsec - begin_grid.tv_nsec) / NSEC_SEC_MUL +
                                            (end_main.tv_sec - begin_grid.tv_sec));

    return EXIT_SUCCESS;
}

// grid search function with loop variables

void gridloopsearch(
    double dd1, double dd2, double dd3, double dd4, double dd5, double dd6, double dd7, double dd8,
    double dd9, double dd10, double dd11, double dd12, double dd13, double dd14, double dd15,
    double dd16, double dd17, double dd18, double dd19, double dd20, double dd21, double dd22,
    double dd23, double dd24, double dd25, double dd26, double dd27, double dd28, double dd29,
    double dd30, double c11, double c12, double c13, double c14, double c15, double c16, double c17,
    double c18, double c19, double c110, double d1, double ey1, double c21, double c22, double c23,
    double c24, double c25, double c26, double c27, double c28, double c29, double c210, double d2,
    double ey2, double c31, double c32, double c33, double c34, double c35, double c36, double c37,
    double c38, double c39, double c310, double d3, double ey3, double c41, double c42, double c43,
    double c44, double c45, double c46, double c47, double c48, double c49, double c410, double d4,
    double ey4, double c51, double c52, double c53, double c54, double c55, double c56, double c57,
    double c58, double c59, double c510, double d5, double ey5, double c61, double c62, double c63,
    double c64, double c65, double c66, double c67, double c68, double c69, double c610, double d6,
    double ey6, double c71, double c72, double c73, double c74, double c75, double c76, double c77,
    double c78, double c79, double c710, double d7, double ey7, double c81, double c82, double c83,
    double c84, double c85, double c86, double c87, double c88, double c89, double c810, double d8,
    double ey8, double c91, double c92, double c93, double c94, double c95, double c96, double c97,
    double c98, double c99, double c910, double d9, double ey9, double c101, double c102,
    double c103, double c104, double c105, double c106, double c107, double c108, double c109,
    double c1010, double d10, double ey10, double kk)
{
    // results values
    double x2, x3, x4, x5, x6, x7, x8, x9, x10;

    // constraint values
    double q1, q2, q3, q4, q5, q6, q7, q8, q9, q10;

    // results points
    long pnts = 0;

    // re-calculated limits

    // opening the "results-v2.txt" for writing he results in append mode
    FILE *fptr = fopen("./results-v2.txt", "w");
    if (fptr == NULL)
    {
        printf("Error in creating file !");
        exit(1);
    }

    // initialization of re calculated limits, xi's.
    const double e1 = kk * ey1;
    const double e2 = kk * ey2;
    const double e3 = kk * ey3;
    const double e4 = kk * ey4;
    const double e5 = kk * ey5;
    const double e6 = kk * ey6;
    const double e7 = kk * ey7;
    const double e8 = kk * ey8;
    const double e9 = kk * ey9;
    const double e10 = kk * ey10;

    x2 = dd4 - dd6;
    x3 = dd7 - dd9;
    x4 = dd10 - dd12;
    x5 = dd13 - dd15;
    x6 = dd16 - dd18;
    x7 = dd19 - dd21;
    x8 = dd22 - dd24;
    x9 = dd25 - dd27;
    x10 = dd28 - dd30;

    // for loop upper values
    const int s1 = floor((dd2 - dd1) / dd3);
    const int s2 = floor((dd5 - dd4) / dd6);
    const int s3 = floor((dd8 - dd7) / dd9);
    const int s4 = floor((dd11 - dd10) / dd12);
    const int s5 = floor((dd14 - dd13) / dd15);
    const int s6 = floor((dd17 - dd16) / dd18);
    const int s7 = floor((dd20 - dd19) / dd21);
    const int s8 = floor((dd23 - dd22) / dd24);
    const int s9 = floor((dd26 - dd25) / dd27);
    const int s10 = floor((dd29 - dd28) / dd30);

    // grid search starts
    
#pragma omp parallel for
    for (int r1 = 0; r1 < s1; ++r1)
    {
        double x1 = dd1 + r1 * dd3;
        double x2 = dd4 - dd6;
        double q1_1 = c11 * x1 - d1;
        double q2_1 = c21 * x1 - d2;
        double q3_1 = c31 * x1 - d3;
        double q4_1 = c41 * x1 - d4;
        double q5_1 = c51 * x1 - d5;
        double q6_1 = c61 * x1 - d6;
        double q7_1 = c71 * x1 - d7;
        double q8_1 = c81 * x1 - d8;
        double q9_1 = c91 * x1 - d9;
        double q10_1 = c101 * x1 - d10;

        for (int r2 = 0; r2 < s2; ++r2)
        {
            x2 += dd6;
            double x3 = dd7 - dd9;
            double q1_2 = q1_1 + c12 * x2;
            double q2_2 = q2_1 + c22 * x2;
            double q3_2 = q3_1 + c32 * x2;
            double q4_2 = q4_1 + c42 * x2;
            double q5_2 = q5_1 + c52 * x2;
            double q6_2 = q6_1 + c62 * x2;
            double q7_2 = q7_1 + c72 * x2;
            double q8_2 = q8_1 + c82 * x2;
            double q9_2 = q9_1 + c92 * x2;
            double q10_2 = q10_1 + c102 * x2;

            for (int r3 = 0; r3 < s3; ++r3)
            {
                x3 += dd9;
                double q1_3 = q1_2 + c13 * x3;
                double q2_3 = q2_2 + c23 * x3;
                double q3_3 = q3_2 + c33 * x3;
                double q4_3 = q4_2 + c43 * x3;
                double q5_3 = q5_2 + c53 * x3;
                double q6_3 = q6_2 + c63 * x3;
                double q7_3 = q7_2 + c73 * x3;
                double q8_3 = q8_2 + c83 * x3;
                double q9_3 = q9_2 + c93 * x3;
                double q10_3 = q10_2 + c103 * x3;

                for (int r4 = 0; r4 < s4; ++r4)
                {
                    x4 = dd10 + r4 * dd12;
                    double x5 = dd13 - dd15;
                    double q1_4 = q1_3 + c14 * x4;
                    double q2_4 = q2_3 + c24 * x4;
                    double q3_4 = q3_3 + c34 * x4;
                    double q4_4 = q4_3 + c44 * x4;
                    double q5_4 = q5_3 + c54 * x4;
                    double q6_4 = q6_3 + c64 * x4;
                    double q7_4 = q7_3 + c74 * x4;
                    double q8_4 = q8_3 + c84 * x4;
                    double q9_4 = q9_3 + c94 * x4;
                    double q10_4 = q10_3 + c104 * x4;

                    for (int r5 = 0; r5 < s5; ++r5)
                    {
                        x5 += dd15;
                        double x6 = dd16 - dd18;
                        double q1_5 = q1_4 + c15 * x5;
                        double q2_5 = q2_4 + c25 * x5;
                        double q3_5 = q3_4 + c35 * x5;
                        double q4_5 = q4_4 + c45 * x5;
                        double q5_5 = q5_4 + c55 * x5;
                        double q6_5 = q6_4 + c65 * x5;
                        double q7_5 = q7_4 + c75 * x5;
                        double q8_5 = q8_4 + c85 * x5;
                        double q9_5 = q9_4 + c95 * x5;
                        double q10_5 = q10_4 + c105 * x5;

                        for (int r6 = 0; r6 < s6; ++r6)
                        {
                            x6 += dd18;
                            double x7 = dd19 - dd21;
                            double q1_6 = q1_5 + c16 * x6;
                            double q2_6 = q2_5 + c26 * x6;
                            double q3_6 = q3_5 + c36 * x6;
                            double q4_6 = q4_5 + c46 * x6;
                            double q5_6 = q5_5 + c56 * x6;
                            double q6_6 = q6_5 + c66 * x6;
                            double q7_6 = q7_5 + c76 * x6;
                            double q8_6 = q8_5 + c86 * x6;
                            double q9_6 = q9_5 + c96 * x6;
                            double q10_6 = q10_5 + c106 * x6;

                            for (int r7 = 0; r7 < s7; ++r7)
                            {
                                x7 += dd21;
                                double x8 = dd22 - dd24;
                                double q1_7 = q1_6 + c17 * x7;
                                double q2_7 = q2_6 + c27 * x7;
                                double q3_7 = q3_6 + c37 * x7;
                                double q4_7 = q4_6 + c47 * x7;
                                double q5_7 = q5_6 + c57 * x7;
                                double q6_7 = q6_6 + c67 * x7;
                                double q7_7 = q7_6 + c77 * x7;
                                double q8_7 = q8_6 + c87 * x7;
                                double q9_7 = q9_6 + c97 * x7;
                                double q10_7 = q10_6 + c107 * x7;
                                for (int r8 = 0; r8 < s8; ++r8)
                                {
                                    x8 += dd24;
                                    double x9 = dd25 - dd27;
                                    double q1_8 = q1_7 + c18 * x8;
                                    double q2_8 = q2_7 + c28 * x8;
                                    double q3_8 = q3_7 + c38 * x8;
                                    double q4_8 = q4_7 + c48 * x8;
                                    double q5_8 = q5_7 + c58 * x8;
                                    double q6_8 = q6_7 + c68 * x8;
                                    double q7_8 = q7_7 + c78 * x8;
                                    double q8_8 = q8_7 + c88 * x8;
                                    double q9_8 = q9_7 + c98 * x8;
                                    double q10_8 = q10_7 + c108 * x8;

                                    for (int r9 = 0; r9 < s9; ++r9)
                                    {
                                        x9 += dd27;
                                        double x10 = dd28 - dd30;
                                        double q1_9 = q1_8 + c19 * x9;
                                        double q2_9 = q2_8 + c29 * x9;
                                        double q3_9 = q3_8 + c39 * x9;
                                        double q4_9 = q4_8 + c49 * x9;
                                        double q5_9 = q5_8 + c59 * x9;
                                        double q6_9 = q6_8 + c69 * x9;
                                        double q7_9 = q7_8 + c79 * x9;
                                        double q8_9 = q8_8 + c89 * x9;
                                        double q9_9 = q9_8 + c99 * x9;
                                        double q10_9 = q10_8 + c109 * x9;

                                        for (int r10 = 0; r10 < s10; ++r10)
                                        {
                                            x10 += dd30;

                                            // constraints

                                            if ((fabs(q1_9 + c110 * x10) <= e1) && (fabs(q2_9 + c210 * x10) <= e2) && (fabs(q3_9 + c310 * x10) <= e3) && (fabs(q4_9 + c410 * x10) <= e4) && (fabs(q5_9 + c510 * x10) <= e5) &&
                                                (fabs(q6_9 + c610 * x10) <= e6) && (fabs(q7_9 + c710 * x10) <= e7) && (fabs(q8_9 + c810 * x10) <= e8) && (fabs(q9_9 + c910 * x10) <= e9) && (fabs(q10_9 + c1010 * x10) <= e10))
                                            {

#pragma omp critical
                                                {
                                                    pnts = pnts + 1;
                                                    // xi's which satisfy the constraints to be written in file
                                                    fprintf(fptr, "%lf\t", x1);
                                                    fprintf(fptr, "%lf\t", x2);
                                                    fprintf(fptr, "%lf\t", x3);
                                                    fprintf(fptr, "%lf\t", x4);
                                                    fprintf(fptr, "%lf\t", x5);
                                                    fprintf(fptr, "%lf\t", x6);
                                                   fprintf(fptr, "%lf\t", x7);
                                                    fprintf(fptr, "%lf\t", x8);
                                                    fprintf(fptr, "%lf\t", x9);
                                                    fprintf(fptr, "%lf\n", x10);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    fclose(fptr);
    printf("result pnts: %ld\n", pnts);

    // end function gridloopsearch
}

