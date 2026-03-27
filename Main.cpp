#include <iostream>
#include <fstream>
#include <string>
#include <cmath> 

using namespace std;

struct Node {
    double x; 
    double y; 
    bool BC = false; 
};

struct Element {
	int ID[4]; 
};

struct Jakobian {
	double J[2][2];
	double detJ; 
    double invJ[2][2]; 
};


struct ElemUniv {
	double* dn_dksi; 
    double* dn_deta; 
};

struct Grid {
    int numNodes; 
	int numElements; 
    Node* nodes;  
	Element* elements; 
};

struct GlobalData {
    double simulationTime;
    double timeStep; 
	double conductivity; 
    double alfa; 
	double tot; 
	double initialTemperature; 
    double density; 
	double specificHeat; 
	int numNodes; 
	int numElements; 
	int npc; 
};


static void get2P(double* n, double* w) {
    double val = 1.0 / sqrt(3.0);
    n[0] = -val;
    n[1] = val;
    w[0] = 1.0;
    w[1] = 1.0;
}

static void get3P(double* n, double* w) {
    double val1 = sqrt(3.0 / 5.0);
    n[0] = -val1;
    n[1] = 0.0;
    n[2] = val1;
    w[0] = 5.0 / 9.0;
    w[1] = 8.0 / 9.0;
    w[2] = 5.0 / 9.0;
}

static void get4P(double* n, double* w) {
    double val1 = 0.861136;
    double val2 = 0.339981;
    n[0] = -val1;
    n[1] = -val2;
    n[2] = val2;
    n[3] = val1;
    double weight1 = 0.347855;
    double weight2 = 0.652145;
    w[0] = weight1;
    w[1] = weight2;
    w[2] = weight2;
    w[3] = weight1;
}


static void funkcjeKsztaltu(double ksi, double eta, double* N) {
    N[0] = 0.25 * (1 - ksi) * (1 - eta);
    N[1] = 0.25 * (1 + ksi) * (1 - eta);
    N[2] = 0.25 * (1 + ksi) * (1 + eta);
    N[3] = 0.25 * (1 - ksi) * (1 + eta);
}


static void pochodne(double ksi, double eta, double* dN_dKsi, double* dN_dEta) {
    dN_dKsi[0] = -0.25 * (1 - eta);  
    dN_dKsi[1] = 0.25 * (1 - eta);   
    dN_dKsi[2] = 0.25 * (1 + eta);   
    dN_dKsi[3] = -0.25 * (1 + eta);  

    dN_dEta[0] = -0.25 * (1 - ksi);  
    dN_dEta[1] = -0.25 * (1 + ksi);  
    dN_dEta[2] = 0.25 * (1 + ksi);   
    dN_dEta[3] = 0.25 * (1 - ksi);   
}


void solveGauss(double** A, double* B, double* X, int n) { 

	
    for (int i = 0; i < n; i++) {
		int max = i; 
        for (int k = i + 1; k < n; k++)
            if (abs(A[k][i]) > abs(A[max][i])) max = k;

        
        double* tempA = A[i]; A[i] = A[max]; A[max] = tempA;
        double tempB = B[i]; B[i] = B[max]; B[max] = tempB;
		
        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i] / A[i][i];
            B[k] -= factor * B[i];
            for (int j = i; j < n; j++) A[k][j] -= factor * A[i][j];
        }
    }
   
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) sum += A[i][j] * X[j];
        X[i] = (B[i] - sum) / A[i][i];
    }
}


int main() {
    
    GlobalData gd;
    Grid grid;

    ifstream file("dane2.txt");
    if (!file.is_open()) {
        cout << "Nie udalo sie otworzyc pliku!" << endl;
        return 1;
    }

    string temp;
    file >> temp >> gd.simulationTime;
    file >> temp >> gd.timeStep;
    file >> temp >> gd.conductivity;
    file >> temp >> gd.alfa;
    file >> temp >> gd.tot;
    file >> temp >> gd.initialTemperature;
    file >> temp >> gd.density;
    file >> temp >> gd.specificHeat;
    file >> temp >> temp >> gd.numNodes;
    file >> temp >> temp >> gd.numElements;
	//file >> temp >> gd.npc; 
	gd.npc = 4; 

    grid.numNodes = gd.numNodes;
    grid.numElements = gd.numElements;
    grid.nodes = new Node[grid.numNodes];
    grid.elements = new Element[grid.numElements];

    file >> temp; 
    for (int i = 0; i < grid.numNodes; i++) {
        int id;
        char comma;
        file >> id >> comma >> grid.nodes[i].x >> comma >> grid.nodes[i].y;
    }

    file >> temp >> temp; 
    for (int i = 0; i < grid.numElements; i++) {
        int id;
        char comma;
        file >> id >> comma
            >> grid.elements[i].ID[0] >> comma
            >> grid.elements[i].ID[1] >> comma
            >> grid.elements[i].ID[2] >> comma
            >> grid.elements[i].ID[3];
    }

    file >> temp; 
    int bcID;
    char comma;
    while (file >> bcID) {
        grid.nodes[bcID - 1].BC = true; 
        if (!(file >> comma)) break;    
    }

    file.close();
    
    int nGauss = gd.npc;
    double pGauss[10] = { 0 };
    double wGauss[10] = { 0 };

    if (nGauss == 2) get2P(pGauss, wGauss);
    else if (nGauss == 3) get3P(pGauss, wGauss);
    else if (nGauss == 4) get4P(pGauss, wGauss);

    
    double** H = new double* [grid.numNodes];
    for (int i = 0; i < grid.numNodes; i++) {
        H[i] = new double[grid.numNodes];
        for (int j = 0; j < grid.numNodes; j++) H[i][j] = 0.0; 
    }

    double* P = new double[grid.numNodes];
    for (int i = 0; i < grid.numNodes; i++) P[i] = 0.0;

    
    double** C = new double* [grid.numNodes];
    for (int i = 0; i < grid.numNodes; i++) {
        C[i] = new double[grid.numNodes];
        for (int j = 0; j < grid.numNodes; j++) C[i][j] = 0.0;
    }
    
    
    for (int e = 0; e < grid.numElements; e++) {
        //cout << "Element " << e + 1 << ":" << endl;

        Jakobian jak;
        double Hlocal[4][4] = { 0 };
        double Hbc[4][4] = { 0 };
        double Plocal[4] = { 0 };
        double Clocal[4][4] = { 0 };

        
        for (int i = 0; i < nGauss; i++) {
            for (int j = 0; j < nGauss; j++) {

                double dN_dksi[4], dN_deta[4];
                pochodne(pGauss[i], pGauss[j], dN_dksi, dN_deta);
                double N[4];
                funkcjeKsztaltu(pGauss[i], pGauss[j], N);


                jak.J[0][0] = 0;
                jak.J[0][1] = 0;
                jak.J[1][0] = 0;
                jak.J[1][1] = 0;

                for (int k = 0; k < 4; k++) {
                    int nodeIdx = grid.elements[e].ID[k] - 1; 
                    double nodeX = grid.nodes[nodeIdx].x;
                    double nodeY = grid.nodes[nodeIdx].y;

                    jak.J[0][0] += dN_dksi[k] * nodeX; 
                    jak.J[0][1] += dN_dksi[k] * nodeY; 
                    jak.J[1][0] += dN_deta[k] * nodeX; 
                    jak.J[1][1] += dN_deta[k] * nodeY; 
                }

                
                jak.detJ = jak.J[0][0] * jak.J[1][1] - jak.J[0][1] * jak.J[1][0];

                
                jak.invJ[0][0] = jak.J[1][1] / jak.detJ;
                jak.invJ[0][1] = -jak.J[0][1] / jak.detJ;
                jak.invJ[1][0] = -jak.J[1][0] / jak.detJ;
                jak.invJ[1][1] = jak.J[0][0] / jak.detJ;

                /*cout << "  Punkt PC(" << i + 1 << "," << j + 1 << "):" << endl;

               
                cout << "    Macierz J:    [" << jak.J[0][0] << " " << jak.J[0][1] << "]" << endl;
                cout << "                  [" << jak.J[1][0] << " " << jak.J[1][1] << "]" << endl;

                
                cout << "    Wyznacznik:   " << jak.detJ << endl;

                
                cout << "    Macierz invJ: [" << jak.invJ[0][0] << " " << jak.invJ[0][1] << "]" << endl;
                cout << "                  [" << jak.invJ[1][0] << " " << jak.invJ[1][1] << "]" << endl << endl;*/


                double dN_dx[4], dN_dy[4];
                for (int k = 0; k < 4; k++) {
                    dN_dx[k] = jak.invJ[0][0] * dN_dksi[k] + jak.invJ[0][1] * dN_deta[k];
                    dN_dy[k] = jak.invJ[1][0] * dN_dksi[k] + jak.invJ[1][1] * dN_deta[k];
                }

                double weight = wGauss[i] * wGauss[j];

     
                for (int row = 0; row < 4; row++) {
                    for (int col = 0; col < 4; col++) {
                        double temp = (dN_dx[row] * dN_dx[col] + dN_dy[row] * dN_dy[col]);
                        Hlocal[row][col] += gd.conductivity * temp * jak.detJ * weight;
                    }
                }

                for (int r = 0; r < 4; r++) {
                    for (int c = 0; c < 4; c++) {
                        Clocal[r][c] += gd.density * gd.specificHeat * N[r] * N[c] * jak.detJ * weight;
                    }
                }

            }
        }

        /*cout << "  Macierz H lokalna dla elementu:" << e + 1 << endl;
        for (int row = 0; row < 4; row++) {
            cout << "    ";
            for (int col = 0; col < 4; col++) {
                cout << Hlocal[row][col] << "\t";
            }
            cout << endl;
        }
        cout << endl;*/

        int boki[4][2] = { {0,1}, {1,2}, {2,3}, {3,0} };

        for (int b = 0; b < 4; b++) {
            int node1_loc = boki[b][0];
            int node2_loc = boki[b][1];
            int idx1 = grid.elements[e].ID[node1_loc] - 1;
            int idx2 = grid.elements[e].ID[node2_loc] - 1;


            if (grid.nodes[idx1].BC && grid.nodes[idx2].BC) {
    
                double L = sqrt(pow(grid.nodes[idx2].x - grid.nodes[idx1].x, 2) + pow(grid.nodes[idx2].y - grid.nodes[idx1].y, 2));
                double detJ_bc = L / 2.0;

                for (int p = 0; p < nGauss; p++) {
                    double N[4];
                    double ksi, eta;

                    if (b == 0) { ksi = pGauss[p]; eta = -1.0; }
                    else if (b == 1) { ksi = 1.0; eta = pGauss[p]; }
                    else if (b == 2) { ksi = -pGauss[p]; eta = 1.0; }
                    else { ksi = -1.0; eta = -pGauss[p]; }

                    funkcjeKsztaltu(ksi, eta, N);

                    for (int row = 0; row < 4; row++) {
                        for (int col = 0; col < 4; col++) {
       
                            Hbc[row][col] += gd.alfa * N[row] * N[col] * wGauss[p] * detJ_bc;
                        }
                    }

 
                    for (int i_local = 0; i_local < 4; i_local++) {
                        Plocal[i_local] += gd.alfa * N[i_local] * gd.tot * wGauss[p] * detJ_bc;
                    }

                    
                }

            }
            /*cout << "  Macierz Hbc lokalna dla elementu:" << e + 1 << endl;
            for (int row = 0; row < 4; row++) {
                cout << "    ";
                for (int col = 0; col < 4; col++) {
                    cout << Hbc[row][col] << "\t";
                }
                cout << endl;
            }
            cout << endl;

            cout << "  Wektor P lokalny dla elementu:" << e + 1 << endl;
            for (int row = 0; row < 4; row++) {
                cout << "    ";
                cout << Plocal[row] << "\t";
                cout << endl;
            }
            cout << endl;

            cout << "  Macierz C lokalna dla elementu:" << e + 1 << endl;
            for (int row = 0; row < 4; row++) {
                cout << "    ";
                for (int col = 0; col < 4; col++) {
                    cout << Clocal[row][col] << "\t";
                }
                cout << endl;
            }
            cout << endl;*/
        }


        for (int i = 0; i < 4; i++) {
            int rG = grid.elements[e].ID[i] - 1;
            P[rG] += Plocal[i];
            for (int j = 0; j < 4; j++) {
                int cG = grid.elements[e].ID[j] - 1;
                H[rG][cG] += (Hlocal[i][j] + Hbc[i][j]);
                C[rG][cG] += Clocal[i][j];
            }
        }
    }




    /*cout << "MACIERZ GLOBALNA H:" << endl;
    for (int i = 0; i < grid.numNodes; i++) {
        for (int j = 0; j < grid.numNodes; j++) {
            cout << H[i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl;

    cout << "WEKTOR GLOBALNY P:" << endl;
    for (int i = 0; i < grid.numNodes; i++) {
        cout << P[i] << "\t";
        cout << endl;
    }
    cout << endl;

    cout << "MACIERZ GLOBALNA C:" << endl;
    for (int i = 0; i < grid.numNodes; i++) {
        for (int j = 0; j < grid.numNodes; j++) {
            cout << C[i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl;


    cout << "WSPOLRZEDNE WEZLOW:" << endl;
    for (int i = 0; i < grid.numNodes; i++) {
        cout << "Node " << i + 1 << ": x=" << grid.nodes[i].x << " y=" << grid.nodes[i].y << endl;
    }

    cout << "\nELEMENTY (ID WEZLOW):" << endl;
    for (int i = 0; i < grid.numElements; i++) {
        cout << "Element " << i + 1 << ": ";
        for (int j = 0; j < 4; j++) {
            cout << grid.elements[i].ID[j] << " ";
        }
        cout << endl;
    }

    cout << "\nWLASCIWOSCI GLOBALNE:" << endl;
    cout << "Czas symulacji: " << gd.simulationTime << endl;
    cout << "Krok czasowy: " << gd.timeStep << endl;
    cout << "Przewodnosc: " << gd.conductivity << endl;
    cout << "Wspolczynnik alfa: " << gd.alfa << endl;
    cout << "Temperatura otoczenia: " << gd.tot << endl;
    cout << "Poczatkowa temperatura: " << gd.initialTemperature << endl;
    cout << "Gestosc: " << gd.density << endl;
    cout << "Cieplo wlasciwe: " << gd.specificHeat << endl;
    cout << endl;*/


    double* t0 = new double[grid.numNodes]; 
	double* t1 = new double[grid.numNodes]; 
	
    for (int i = 0; i < grid.numNodes; i++) t0[i] = gd.initialTemperature;

    double dT = gd.timeStep;
	for (double time = dT; time <= gd.simulationTime; time += dT) { 
        
        double** A = new double* [grid.numNodes];
        double* B = new double[grid.numNodes];
        for (int i = 0; i < grid.numNodes; i++) {
            A[i] = new double[grid.numNodes];
            B[i] = 0;
            for (int j = 0; j < grid.numNodes; j++) {
                A[i][j] = H[i][j] + (C[i][j] / dT); 
                B[i] += (C[i][j] / dT) * t0[j];
            }
            B[i] += P[i];
        }

        solveGauss(A, B, t1, grid.numNodes);

        double minT = t1[0], maxT = t1[0];
        for (int i = 0; i < grid.numNodes; i++) {
            if (t1[i] < minT) minT = t1[i];
            if (t1[i] > maxT) maxT = t1[i];
            t0[i] = t1[i]; 
        }
        cout << "Time: " << time << "s \t Min: " << minT << " \t Max: " << maxT << endl;

        
        for (int i = 0; i < grid.numNodes; i++) delete[] A[i];
        delete[] A; delete[] B;
    }

    


    for (int i = 0; i < grid.numNodes; i++) {
        delete[] H[i];
        delete[] C[i];
    }
	delete[] H;
	delete[] P;
	delete[] t0;
	delete[] t1;   
	delete[] C;
    delete[] grid.nodes;
    delete[] grid.elements;

    return 0;
}