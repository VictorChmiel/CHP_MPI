#include <iostream>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <time.h>

/*  Ecole des Mines de Saint-Étienne - avril 2018 *
 *  @ J. Bruchon, N. Moulin
 *  Toolbox CHP - UP MPI
 * 
 * Résolution de Laplacien(u) = 0 par un algorithme de Monte-Carlo */

using namespace std;


int main(int argc, char ** argv)
{
	int rang, nb_procs;

	// Lecture des données d'entrée dans le fichier passé en paramètre
	

	if ( argc < 2 )
	{
		if ( rang== 0 )
		  cout << "Il manque le nom du fichier de données !" << endl;
		exit(1);
        }

	int nx, ny, nb_tirages;
	double temps_init, temps_final;

	// Lecture du fichier de données
	

	ifstream fic_donnees;
	fic_donnees.open(argv[1], ifstream::in);

        fic_donnees >> nx >> ny >> nb_tirages;
	fic_donnees.close();

	// Affichage des données lues : permet de vérifier que tout va bien.

	if ( rang == 0 )
	{
		cout << "************************" << endl;
		cout << "Calcul exécuté sur " << nb_procs << " coeurs" << endl;
		cout << "Dimension selon x : " << nx << endl;
		cout << "Dimension selon y : " << ny << endl;
		cout << "Nombre de lancés par case : " << nb_tirages << endl;
	}

	// Initialisation et allocation
	
	double * grille = new double [nx*ny];
	int i,j;

	double valeurs_aux_bords[4] = {0,1,2,3}; // valeurs de la solution sur les bords y=0, x=nx, y=ny et x=0
	
	for (i=0; i<nx*ny; i++)
		grille[i] = 0;  // Initialisation de la solution
			
	// Initialisation de la solution sur les bords
	
	for (i=0; i<nx; i++) // y=0
	  grille[i] = valeurs_aux_bords[0];
	
	for (i=0; i<ny; i++) // x=nx
	  grille[i*nx + nx-1] = valeurs_aux_bords[1];
	
	for (i=0; i<nx; i++) // y=ny
	  grille[(nx-1)*ny + i] = valeurs_aux_bords[2];

	for (i=0; i<ny; i++) // x=0
	  grille[i*nx] = valeurs_aux_bords[3];

	srand(time(NULL)); // Initialisation du générateur de nombres aléatoire

	// 1. Boucles sur j et i
	
	for (j=1; j<ny-1; j++)
	{
	  for (i=1; i<nx-1; i++) 
	  {
	    for (int n=0; n<nb_tirages; n++) // on effectue nb_tirages tirages de particules
	    {
	      int position[2] = {i,j}; // Initialisation de la position de la particule
	      bool stop = 0; // variable booléenne qui passe à 1 quand on rencontre un bord
	      double valeur = 0;
	      while ( stop == 0 ) // tant que stop est nul, faire ...
	      {
		  int decision = rand() % 4;  // nombre aléatoire [0,1]
		  
		  switch(decision){

			case 0:
				position[0]--;
			break;

		  	case 1: 
				position[1]++;
			break;

		      	case 2: 
				position[0]++;
			break;

			case 3: 
				position[1]--;
			break;
		}

		  // Test des bords : est-on sur un bord, si oui lequel ?

		  if ( ( position[0] == 0 ) || ( position[0] == nx-1 ) || ( position[1] == 0 ) || ( position[1] == ny-1 ) )
		  {
		    valeur = grille[position[1]*nx+ position[0]]; // On récupère la valeur au bord atteint
		    stop = 1;
		  }
	      }
	      
	      grille[j*nx+i] += valeur;
	    }
	    
	    // Valeur finale de la solution en (i,j) : c'est la moyenne des tirages
	    
	    grille[j*nx+i] /= nb_tirages;
	  }
	}
	
	// 2. Ecriture au format vtk

	cout << "Écriture du fichier vtk" << endl;
	stringstream nom;
	nom << "stochastique";
	int Nbnoe = nx*ny;
	double dx = 1.0/(nx-1);
	double dy = 1.0/(ny-1);

	nom << ".vtk" << '\0';
	ofstream fic(nom.str().c_str());

	fic << "# vtk DataFile Version 2.0" << endl;
	fic << "Laplacien stochastique" << endl;
	fic << "ASCII" << endl;
	fic << "DATASET STRUCTURED_POINTS" << endl;
	fic << "DIMENSIONS " << nx << "  " << ny << "  1 " << endl;
	fic << "ORIGIN 0 0 0" << endl;
	fic << "SPACING " << dx << "  " << dy << "  1" << endl;
	fic << "POINT_DATA " << Nbnoe << endl;
	fic << "SCALARS Concentration float" << endl;
	fic << "LOOKUP_TABLE default" << endl;
	for (i=0; i<Nbnoe; i++)
	{
	  fic << grille[i] << endl;
	}
	fic.close();				


	// Fin

	delete [] grille;
	return 0;
}


