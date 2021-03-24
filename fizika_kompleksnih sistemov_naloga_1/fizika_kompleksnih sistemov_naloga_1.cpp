

#include <iostream>
#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#pragma warning(disable:4996)
FILE* data2, * out, * data1;

int range_width; //sirina parcele konekcija
const int maxSIZE = 5000;
int i, ii, j, jj, N, e[maxSIZE][maxSIZE], k[maxSIZE], n_node[maxSIZE][maxSIZE], arrConnections[maxSIZE];;
int k_ij;
double p;
double c_0, l_0, cc, ll;
int ij;
int j_new, ok;
double vsota, ef, path;
int e1[maxSIZE][maxSIZE], dist[maxSIZE], oce1[maxSIZE], R, u, cent[maxSIZE][maxSIZE];
double xk[maxSIZE], yk[maxSIZE], dij[maxSIZE][maxSIZE], w_k[maxSIZE], delta;
double c[maxSIZE], Cavg;
double rr, k_avg;
int n, kj_tot, m;
double pki, euclid_dist; 

void interactionsCount(void), regularnetwork(void), smallworld(void), SFnetwork(void), clustering(void), 
shortestpath(void), Spatial_ScaleFree_Network(void), clustering_And_Shortest_Path(void), shortest_path_Nrising(void), euclidian_distance(void), spatial_scalefree_delta(void);
int main()
{
    srand(time(NULL));
    // maloga 1
    /* N = 100;
    p = 0.3;
   // regularnetwork(void);
    //smallworld();
    range_width = 2;
	interactionsCount();
    */
    N = 1000;
    delta = -10.0;
    Spatial_ScaleFree_Network();
    // naloga 2
    //clustering_And_Shortest_Path();
    //naloga 3
    //shortest_path_Nrising();
    //naloga 4
    //spatial_scalefree_delta();
    out = fopen("cpl_matric.txt", "w");
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            fprintf(out, "%i %i %i\n", i, j, e[i][j]);//xyz zapis v matriko
        }
    }
    fclose(out);
    return 0;
}

void interactionsCount(void)
{
	//prolazim kroz celu mrezu
	int max_num = 0; // najveci broj koneckija
	int temp; // privremena variabla - holder
	for (i = 0; i < N; i++) {
		temp = 0;
		for (j = 0; j < N; j++) {
			arrConnections[i] += e[i][j];
			temp += e[i][j];
			//printf("%i\n", temp);
		}
        printf("%i %i\n", i, temp);
		if (temp > max_num)  max_num = temp;
	}
	int x; //broj povezanosti u range_width
	out = fopen("distribution.txt", "w");
	for (i = 0; i < max_num + range_width; i += range_width) {
		x = 0;
		for (j = 0; j < N; j++) {
			if ((arrConnections[j] >= i) && (arrConnections[j] < i + range_width)) {
				x++;
			}
		}
        if(x > 0) fprintf(out, "%lf %i\n", double(i) + double(range_width) / 2.0, x);
		
	}

	fclose(out);
}

void regularnetwork(void) {
	for (i = 0; i < N; i++) {//po vrsticah
		for (j = i; j < N; j++) {//po stolpcih(samo nad diagonalo)
			if (j == (i + 1) || (j == (i + 2))) {//levi sosed od itega vozla
				e[i][j] = 1;
				e[j][i] = 1;//dodaj še desnega soseda
			}
			else {
				e[i][j] = 0;//ni povezave
				e[j][i] = 0;
			}
		}
	}
	//dodamo periodične robne pogoje
	e[0][N - 1] = 1;//prvi z zadnjim
	e[N - 1][0] = 1;//zadnji s prvim
	e[0][N - 2] = 1;//periodični robni pogojo za 4 sosede
	e[1][N - 1] = 1;
	e[N - 2][0] = 1;
	e[N - 1][1] = 1;
}

void smallworld(void) {
	regularnetwork();//najprej rabimo regularno mrežo, potem sledi prevezovanje
	for (i = 0; i < N; i++) {//čez vse vozle
		for (j = i; j < N; j++) {//čez naddiagonalne elemente
			if (e[i][j] == 1) {//če sta i in j vozel povezana
				rr = ((double)(rand())) / ((double)(RAND_MAX));//naključno število med 0 in 1

				if (p > rr) {//preveži le z verjetnostjo p
					e[i][j] = 0;//razveži i in j vozel
					e[j][i] = 0;//dvosmerne povezave



					ok = 0;
					while (ok < 1) {//izbiraj dokler ne najdeš primernega novega soseda
						j_new = (rand() % N);//naključno število med 0 in N (brez N)
						if ((i != j_new) && (e[i][j_new] == 0)) {//ne s samim sabo in med njima ne sme obstajat povezava
							e[i][j_new] = 1;
							e[j_new][i] = 1;
							ok = 1;
						}
					}
				}
			}
		}
	}
}

void SFnetwork(void) {
	//double rr;
	//int n, kj_tot, m;
	//double pki;
	for (i = 0; i < N; i++) {//ustvarimo prazno matriko
		for (j = 0; j < N; j++) {
			e[i][j] = 0;
			e[i][j] = 0;
		}
	}
	printf("1 del");
	m = 2;//začnemo z m povezanimi vozli
	for (i = 0; i < m; i++) {
		for (j = (i + 1); j < m; j++) {
			e[i][j] = 1;
			e[j][i] = 1;
		}
	}
	printf("2 del");
	for (n = m; n < N; n++) {//zanka za rast mreže od m do N
		//printf("%i\n", n);
		kj_tot = 0;//vsota vseh povezav
		for (i = 0; i < n; i++) {//za vse dosedanje vozle preštejo število povezav
			k[i] = 0;//število povezav i-tega vozla
			for (j = 0; j < n; j++) {//po i-tem stolpcu v matriki preštejemo sosede
				k[i] = k[i] + e[i][j];
			}
			kj_tot = kj_tot + k[i];//vsota vseh povezav v mreži (od vseh dosedanjih vozlov)
		}
		for (i = 0; i < m; i++) {//vsak novi n-ti vozel povežemo z m dosedanjimi vozli
			rr = ((double)rand() / (double)RAND_MAX);//naključno število med 0 in 1
			//printf("%lf\n", rr);
			pki = 0.0;//verjetnost za povezavo
			j = 0;
			while (pki <= rr) {
				pki = pki + double(k[j]) / double(kj_tot);//sešteva normirane vrednosti za povezavo
				if (pki > rr) {//če je izpoljnjen pogoj za povezovanje
					if (e[n][j] < 1) {//preveri če nista povezana
						e[n][j] = 1;
						e[j][n] = 1;
					}
					else {
						pki = 0.0;//resetiraj
						rr = ((double)rand() / (double)RAND_MAX);//naključno število med 0 in 1
						j = 0;
					}
				}
				j += 1;
			}
		}
	}
}

void Spatial_ScaleFree_Network(void) {
    int II;
    double xkk, ykk, dxy, mindist, w_k_tot;
    int ok1;
    for (i = 0; i < N; i++) {  // ustvari prazno sklopitveno NxN matriko 
        for (j = 0; j < N; j++) {
            e[i][j] = 0;
            e[j][i] = 0;
        }
    }
    // spatial COORDINATES 
    for (int I = 0; I < N; I++) {  // doloci koordinate celic, znotraj enotskega kvadrata 
        xkk = ((double)(rand()) + 1.) / ((double)(RAND_MAX)+1.);
        ykk = ((double)(rand()) + 1.) / ((double)(RAND_MAX)+1.);
        ok1 = 1;
        for (int II = 0; II < I; II++) { 
            // PREDOLGO TRAJA PRI VELIKEM ŠTEVILU IGRALCEV, izkljuciti varovalko za koncne izracune 
            dxy = sqrt((xk[II] - xkk) * (xk[II] - xkk) + (yk[II] - ykk) * (yk[II] - ykk));
            if (dxy < 0.5 * 1.0 / sqrt((double)N)) { 
                // min. razdalja med vozloma, pomembno le za vizualizacijo 
                ok1 = 0;
            }
        }
        if (ok1 > 0) {
            xk[I] = xkk;
            yk[I] = ykk;
        }
        else {
            I = I - 1;
        }
    }
    
    out = fopen("node_coordinates.txt", "w");
    for (i = 0; i < N; i++)
    {
        fprintf(out, "%.4lf %.4lf\n", xk[i], yk[i]);
    }
    fclose(out);
    for (i = 0; i < N; i++)  {
        // izracunaj razdalje med pari celic celic 
        for (j = i; j < N; j++) {
            //razdalje normiramo med 0 in 1
            dij[i][j] = sqrt((xk[i] - xk[j]) * (xk[i] - xk[j]) + (yk[i] - yk[j]) * (yk[i] - yk[j])) / sqrt(2.0);
            dij[j][i] = dij[i][j];
        }
    }
    int m = 2; // zacni z dvema vozloma, !!!!!!   DRUGACE PRI "SPATIAL" NE GRE !!!!!!!!!!!!!!!!!!! 
    for (i = 0; i < m; i++) {
        //  ki sta med sabo povezana 
        for (j = i + 1; j < m; j++) {
            e[i][j] = 1;
            e[j][i] = 1;
        }
    }
    for (int n = m; n < N; n++) {
        // dodajanje vozlov 
        w_k_tot = 0.0;  // skupni degree vseh obstojecih vozlov (zaradi normalizacije) 
        for (i = 0; i < n; i++) {
            // po obstojecih vozlih 
            k[i] = 0;  // degree obstojecih vozlov 
            for (j = 0; j < n; j++) {
                k[i] += e[i][j];
            }
            w_k[i] = k[i] * pow(dij[i][n], delta);
            w_k_tot = w_k_tot + w_k[i];  // weighting with Euclidean distance 
        }
        for (i = 0; i < m; i++) {
            // povezi nov n-ti vozel z m obstojecimi vozli 
            rr = ((double)(rand()) + 1.) / ((double)(RAND_MAX)+1.);
            pki = 0.0;  // verjetnost za povezavo 
            j = 0;
            while (pki <= rr)  // preferencno povezovanje // PREVERI, zamik za !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
            {
                pki = pki + w_k[j] / w_k_tot;  // z deg_tot normiras med 0 in 1 
                if (pki >= rr) {
                    if (e[n][j] == 0) { 
                        // preveri, ce novi vozel n ni ze povezan z IItim 
                        e[n][j] = 1;
                        e[j][n] = 1;
                        int save_j = j;  // memorizes the first connected node to the new node n 
                        w_k_tot = 0.0;
                        for (int i1 = 0; i1 < n; i1++) {
                            k[i1] = 0;
                            for (int j1 = 0; j1 < n; j1++) {
                                k[i1] += e[i1][j1];
                            }
                            w_k[i1] = k[i1] * pow(dij[i1][n], delta);
                            if (i1 == save_j) w_k[i1] = 0.0;
                            w_k_tot = w_k_tot + w_k[i1];
                        }
                    }
                    else  {
                        // ce sta n in II bila povezana ze prej, poskuzi znova 
                        pki = 0.0;
                        rr = ((double)(rand()) + 1.) / ((double)(RAND_MAX)+1.);
                        j = -1;
                    }
                }
                j++;
            }
        }
    }
    /*
    out = fopen("povezave.dat", "w");
    for (i = 0; i < N; i++)  // ustvari prazno sklopitveno NxN matriko
    {
        for (j = 0; j < N; j++)
        {
            if ((i != j) && (e[i][j] > 0))
            {
                fprintf(out, "%.4lf %.4lf %.4lf %.4lf\n", xk[i], yk[i], xk[j], yk[j]);
            }
        }
    }
    fclose(out);*/
}

void shortestpath(void) {
    for (i = 0; i < N; i++)  {
        //shortest path 
        for (j = 0; j < N; j++) {
            if (e[i][j] == 0) {
                e1[i][j] = 10e6;
            }
            else {

                e1[i][j] = e[i][j];
            }
            if (i == j) e1[i][j] = 0;
        }
    }
    vsota = 0.0;
    ef = 0.0;
    for (i = 0; i < N; i++) {
        //shortest path 
       for (j = 0; j < N; j++) {
            dist[j] = e1[i][j];
            oce1[j] = -i - 1;
        }
        dist[i] = 0.0;
        oce1[i] = 0.0;
        for (j = 0; j < N; j++) {
            R = 10e6;
            ok = 1;
            for (ii = 0; ii < N; ii++) {
                if ((oce1[ii] < 0) && (dist[ii] < R)) {
                    R = dist[ii];
                    u = ii;
                }
            }
            if (R == 10e6) {
                ok = 0;
                ii = N + 1;
            }
            if (ok > 0) {
                oce1[u] = -oce1[u];
                for (ii = 0; ii < N; ii++) {
                    if ((oce1[ii] < 0) && (dist[ii] > (dist[u] + e1[u][ii]))) {
                        dist[ii] = dist[u] + e1[u][ii];
                        oce1[ii] = -u - 1;
                    }
                }
            }
        }
        for (j = 0; j < N; j++) {
            vsota = vsota + double(dist[j]);
            cent[i][j] = dist[j];  // shrani najkrajso razdaljo med i-tim in j-tim 
            if (dist[j] == 0) {
                // ce pot ne obstaja 
                ef = ef + 0.0;   // 0 k ucinkovitosti 
            }
            else {
                ef = ef + (1.0 / double(dist[j]));  // vsota recipr. vrednosti najkraj. poti 
            }
        }
    }
    path = vsota / (double(N) * (double(N) - 1.0));
    ef = ef / (double(N) * (double(N) - 1.0));
    /*
    out=fopen("distance_matrix.dat","w");
    for (i=0;i<N;i++) {
        for (j=0;j<N;j++) {
            fprintf(out, "%i %i %i\n", i, j, cent[i][j]);
        }
    }
    fclose(out);*/
}
void clustering(void) {
    //int k_ij;
    for (i = 0; i < N; i++) {
        k[i] = 0;//števec povezav i-tega vozla
        k_ij = 0;//števec števila obstoječih povezav
        for (j = 0; j < N; j++) {//zanka čez vse sosede i-tega vozla
            if (e[i][j] > 0) {
                n_node[i][k[i]] = j;//zapomni si sosede i-tega vozla
                k[i]++;//prešteje št. sosedov i-tega vozla
            }
        }
        for (j = 0; j < k[i]; j++) {//gremo čez vse sosede i-tega vozla
            for (jj = 0; jj < k[i] - j; jj++) {//preveri povezanost vseh j-tih sosedov i-tega vozla
                if (e[n_node[i][j]][n_node[i][j + jj]] > 0) {
                    k_ij++;
                }
            }
        }
        if (k[i] > 1) {
            c[i] = 2.0 * double(k_ij) / (double(k[i]) * (double(k[i]) - 1.0));
        }
        else {
            c[i] = 0.0;
        }
    }
    Cavg = 0.0;//povprečni koeficient gručavosti
    for (i = 0; i < N; i++) {
        Cavg = Cavg + c[i];
    }
    Cavg = Cavg / double(N);



    /*
    out = fopen("clustering.dat", "w");
    for (i = 0; i < N; i++) {
        fprintf(out, "%i %.3lf %.3lf\n", i, c[i], Cavg);
    }
    fclose(out);*/
}

void clustering_And_Shortest_Path(void) {
    out = fopen("shortest_path.txt", "w");
    p = 0.0001;
    N = 100;
    int l = 0;
    double cavg_new = 0.0;
    double path_new = 0.0;
    while (p < 1.0) {
        cavg_new = 0.0;
        path_new = 0.0;
        l = 0;
        while (l < 100) {
            smallworld();
            clustering();
            shortestpath();
            if (path > 2000.0) {
                continue;
            }
            path_new += path;
            cavg_new += Cavg;
            l++;
        }
        path_new = path_new / 100.0;
        cavg_new = cavg_new / 100.0;
        fprintf(out, "%f %f %f\n", p, cavg_new, path_new);
        p *= 1.1;
        printf("%f\n", p);
    }
    fclose(out);
}

void shortest_path_Nrising(void) {
    out = fopen("shortest_path_N.txt", "w");
    p = 0.1;
    N = 100;
    delta = 0.0;
    int l = 0;
    double path_new = 0.0;
    while (N < 1000) {
        path_new = 0.0;
        l = 0;
        while (l < 10) {
            Spatial_ScaleFree_Network();
            shortestpath();
            if (path > 2000.0) {
                continue;
            }
            path_new += path;
            l++;
        }
        path_new = path_new / 10.0;
        fprintf(out, "%i %f\n",N, path_new);
        N *= 1.33;
        printf("%i\n", N);
    }
    fclose(out);
}

void euclidian_distance(void) {
    euclid_dist = 0.0;
    int sum_temp = 0;
    for (i = 0; i < N; i ++) {
        for (j = 0; j < N; j++) {
            if (e[i][j] > 0) {
                euclid_dist += dij[i][j];
                sum_temp++;
            }
        }
    }
    euclid_dist = euclid_dist / double(sum_temp);
}

void spatial_scalefree_delta(void) {
    delta = 0.0;
    N = 500;
    out = fopen("spatial_scalefree_delta.txt", "w");
    while (delta > -10.0) {
        Spatial_ScaleFree_Network();
        clustering();
        shortestpath();
        euclidian_distance();
        if (path > 2000.0) {
            continue;
        }
        fprintf(out, "%lf %lf %lf %lf\n", delta, Cavg, path, euclid_dist);
        delta -= 0.1;
        printf("%lf\n", delta);
    }
    fclose(out);
}