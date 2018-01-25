double F_ST(int ***x,
            int **n,
            int nbr_loci,
            int nbr_demes)

{
    int i, j, k;
    int S_1, S_2;
    int empty_deme;
    double SSI, SSP, MSI, MSP, n_c, n_d;
    double num, den;
    double pi_hat[2];

    num = den = 0.0;
    for (j = 0; j < nbr_loci; ++j)
    {
        S_1 = S_2 = 0;
        empty_deme = 0;
        SSI = SSP = 0.0;
        for (k = 0; k < 2; ++k)
        {
            pi_hat[k] = 0.0;
        }
        for (i = 0; i < nbr_demes; ++i)
        {
            if (n[i][j] == 0)
            {
                empty_deme += 1;
                continue;
            }
            S_1 += n[i][j];
            S_2 += pow(n[i][j], 2);
            for (k = 0; k < 2; ++k)
            {
                pi_hat[k] += (double)x[i][j][k];
                SSI += (double)x[i][j][k] - pow(x[i][j][k], 2) / n[i][j];
            }
        }
        for (k = 0; k < 2; ++k)
        {
            pi_hat[k] /= (double)S_1;
        }
        for (i = 0; i < nbr_demes; ++i)
        {
            if (n[i][j] == 0)
            {
                continue;
            }
            for (k = 0; k < 2; ++k)
            {
                SSP += n[i][j] * pow(((double)x[i][j][k] / n[i][j] -
                                      pi_hat[k]),
                                     2);
            }
        }
        n_d = (nbr_demes - empty_deme);
        n_c = ((double)S_1 - (double)S_2 / S_1) / ((double)n_d - 1.0);
        MSI = SSI / ((double)S_1 - n_d);
        MSP = SSP / ((double)n_d - 1.0);
        num += (MSP - MSI);
        den += (MSP + (n_c - 1.0) * MSI);
    }
    if (den > 0)
    {
        return (num / den);
    }
    else
    {
        return ML_POSINF;
    }
}

double F_ST_poolseq(int ***x,
                    int **n,
                    int nbr_loci,
                    int nbr_demes,
                    int *pool_size)

{
    int i, j, k;
    int R_1, R_2;
    double C_1, C_1_star, n_c;
    double MSI, MSP, SSI, SSP;
    double num, den;
    double pi_hat[2];

    num = den = 0.0;
    for (j = 0; j < nbr_loci; ++j)
    {
        R_1 = R_2 = 0;
        C_1 = C_1_star = 0.0;
        SSI = SSP = 0.0;
        for (k = 0; k < 2; ++k)
        {
            pi_hat[k] = 0.0;
        }
        for (i = 0; i < nbr_demes; ++i)
        {
            C_1 += (double)n[i][j] / pool_size[i] + (double)(pool_size[i] -
                                                             1.0) /
                                                        pool_size[i];
            C_1_star += (double)n[i][j] * ((double)n[i][j] / pool_size[i] +
                                           (double)(pool_size[i] - 1.0) / pool_size[i]);
            R_1 += n[i][j];
            R_2 += n[i][j] * n[i][j];
            if (n[i][j] == 0)
            {
                continue;
            }
            for (k = 0; k < 2; ++k)
            {
                pi_hat[k] += (double)x[i][j][k];
                SSI += (double)x[i][j][k] - pow(x[i][j][k], 2) / n[i][j];
            }
        }
        for (k = 0; k < 2; ++k)
        {
            pi_hat[k] /= (double)R_1;
        }
        C_1_star /= (double)R_1;
        for (i = 0; i < nbr_demes; ++i)
        {
            if (n[i][j] == 0)
            {
                continue;
            }
            for (k = 0; k < 2; ++k)
            {
                SSP += n[i][j] * pow(((double)x[i][j][k] / n[i][j] -
                                      pi_hat[k]),
                                     2);
            }
        }
        n_c = ((double)R_1 - (double)R_2 / R_1) / ((double)C_1 - C_1_star);
        MSI = SSI / ((double)R_1 - C_1);
        MSP = SSP / ((double)C_1 - C_1_star);
        num += (MSP - MSI);
        den += (MSP + (n_c - 1.0) * MSI);
    }
    if (den > 0)
    {
        return (num / den);
    }
    else
    {
        return ML_POSINF;
    }
}
