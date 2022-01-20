#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <math.h>

struct vecteur{
    double *v;
    char *etiq;
    double norme;
};

struct node
{
    double *w;
    char *etiquette;
    double act;
} node;

struct reseau
{
    struct node **r;
} reseau;

typedef struct Bmu{
    int l;
    int c;
    struct Bmu *bmusuiv;
}Bmu;

struct listeBmu{
    struct Bmu *premier;
    int nb;
};

struct vecteur *vecteur;
double *vecMoyen, * vecMin, *vecMax;
int *tabShuffle;
struct listeBmu *liste1;

void init_vecteur(){

    //Allocation de mémoire
    vecteur = malloc(150*sizeof(struct vecteur));
    for(int i=0; i<150; i++){
        vecteur[i].v = malloc(4*sizeof(double));
        vecteur[i].etiq = malloc(30*sizeof(char));
    }

    FILE *fichier;
    fichier = fopen("iris.data","r");
    char chaine[1000];
    char *str, *endPtr;
    const char *delimiter=",\n";
    int i=0, j=0;
    double value, res;

    while (fgets(chaine,1000,fichier)!=NULL){
        str=strtok(chaine,delimiter);

        while(str!=NULL){
            while(j<4){
                value = strtod(str,&endPtr);
                str=strtok(NULL,delimiter);

                vecteur[i].v[j]=value;
                j++;
            }
            j=0;
            // vecteur[i].etiq=str;
            strcpy(vecteur[i].etiq,str);
            str=strtok(NULL,delimiter);

        }
        
        for(int k=0; k<4; k++){
            res+=pow(vecteur[i].v[k],2);
        }
        res = sqrt(res);
        vecteur[i].norme = res;
        res=0.0;

        for(int p=0; p<4; p++){
            vecteur[i].v[p]/=vecteur[i].norme;
            //printf("%.2lf et norm :%.2lf\n",vecteur[i].v[p],vecteur[i].norme);
        }
        i++;
    }
    fclose(fichier);
}

void vecteurMoyen(){
    vecMoyen=malloc(4*sizeof(double));
    double tmpColonne=0.0;
    for(int j=0; j<4; j++){
        for(int i=0; i<150; i++){
            tmpColonne+=vecteur[i].v[j];
        }
        
        vecMoyen[j]= tmpColonne/150;
        //printf("%.2lf\n",vecMoyen[j]);
        tmpColonne=0;
    }
}

void vecteurMinimum(double val){
    vecMin=malloc(4*sizeof(double));
    for(int i=0; i<4; i++){
        vecMin[i]=vecMoyen[i] - val;
        // printf("%.3lf\n",vecMin[i]);
    }
}

void vecteurMaximum(double val){
    vecMax=malloc(4*sizeof(double));
    for(int i=0; i<4; i++){
        vecMax[i]=vecMoyen[i] + val;
        // printf("%.3lf\n",vecMax[i]);
    }
}

void initialisation_shuffle(){
    tabShuffle=malloc(150*sizeof(int));
    for(int i=0; i<150; i++){
        tabShuffle[i]=i;
    }
}

void shuffle(){
    int tmp, alea;
    srand(time(NULL));
    for(int i=0; i<150; i++){
        alea = rand()%150;
        tmp = tabShuffle[i];
        tabShuffle[i]=tabShuffle[alea];
        tabShuffle[alea]=tmp;
    }
}

double *aleatoire(){
    double *vecRes=malloc(4*sizeof(double));
    for(int i=0;i<4;i++){
        vecRes[i]=(vecMax[i]-vecMin[i])*(rand()/(double)RAND_MAX)+vecMin[i];
        // printf("%.2lf\n",vecRes[i]);
    }
    return vecRes;
}

void init_reseauMap(){
    //Allocation de mémoire
    reseau.r=malloc(6*sizeof(struct node*));
    for(int i=0; i<6; i++){
        reseau.r[i]=malloc(10*sizeof(struct node));
    }
    for(int i=0;i<6;i++){
        for(int j=0; j<10;j++){
            reseau.r[i][j].w=malloc(4*sizeof(double));
            reseau.r[i][j].w=aleatoire();
            reseau.r[i][j].etiquette=malloc(30*sizeof(char));
            reseau.r[i][j].etiquette="*";
        }
    }
}

void affichage_reseau(){
    for(int i=0; i<6; i++){
        for(int j=0; j<10; j++){
            printf("%s ", reseau.r[i][j].etiquette);
        }
        printf("\n");
    }
}

struct listeBmu *allocation(){
    struct listeBmu *liste=malloc(sizeof(*liste));
    struct Bmu *bmu = malloc(sizeof(*bmu));
    if(liste==NULL || bmu==NULL){
        exit(EXIT_FAILURE);
    }
    bmu->l=0;
    bmu->c=0;
    bmu->bmusuiv=NULL;
    liste->premier = bmu;
    liste->nb+=1;

    return liste;
}

double distEuclidienne(double *x, double*y, int taille){
    double res=0.0, temp=0.0;
    for(int i=0; i<taille;i++){
        temp =x[i] - y[i];
        res+=pow(temp,2);
    }
    res = sqrt(res);
    return res;
}

void supprimer(struct listeBmu *liste){
    while(liste->premier->bmusuiv){
        struct Bmu *bm = liste->premier;
        liste->premier = liste->premier->bmusuiv;
        free(bm);
    }
    liste->premier=NULL;
    liste->nb =0;
}

void inserer(struct listeBmu *liste, int ligne, int col){
    struct Bmu *bm = malloc(sizeof(Bmu));
    if(liste == NULL || bm ==NULL){ exit(EXIT_FAILURE);}
    bm->l = ligne;
    bm->c=col;
    bm->bmusuiv= liste->premier;
    liste->premier = bm;
    liste->nb = liste->nb + 1;
}

void voisinage(struct Bmu *bmu, double alpha,int ligne, int colonne, int position){
    int rayon=3;
    int x_min, x_max, y_min, y_max;
    for(;rayon>0; rayon--){
        if(bmu->l -rayon<0){
            x_min=0;
        }else{
            x_min= bmu->l - rayon;
        }

        if(bmu->c - rayon <0){
            y_min= 0;
        }else{
            y_min=bmu->c - rayon;
        }

        if(bmu->l+rayon >ligne -1){
            x_max = ligne -1;
        }else{
            x_max=bmu->l + rayon;
        }

        if(bmu->c+rayon >colonne-1){
            y_max=colonne-1;
        }else{
            y_max=bmu->c+rayon;
        }
        for(int i=x_min; i<=x_max;i++){
            for(int j=y_min; j<=y_max; j++){
                for(int k=0; k<4; k++){
                    reseau.r[i][j].w[k]+= (alpha*(vecteur[position].v[k] - reseau.r[i][j].w[k]));
                }
            }
        }
    }
}

void affecte_val(int position,struct Bmu *bm ){
    if(strcmp(vecteur[position].etiq, "Iris-setosa")==0){
        reseau.r[bm->l][bm->c].etiquette="A";
    }else if(strcmp(vecteur[position].etiq, "Iris-versicolor")==0){
        reseau.r[bm->l][bm->c].etiquette="B";
    }else if(strcmp(vecteur[position].etiq, "Iris-virginica")==0){
        reseau.r[bm->l][bm->c].etiquette="C";
    }
}

void apprentissage(int phase, int ligne, int colonne, double alpha){
    double dist,tmp;
    int i,j;
    initialisation_shuffle();
    for(i=0; i<phase; i++){
        shuffle();
        for(j=0; j<150; j++){
            liste1=allocation();
            dist=10.0;
            for(int a=0; a<ligne; a++){
                for(int b=0; b<colonne; b++){
                    tmp=distEuclidienne(vecteur[tabShuffle[j]].v,reseau.r[a][b].w,4);
                    if(tmp<dist){
                        dist=tmp;
                        supprimer(liste1);
                        inserer(liste1,a,b);
                    }
                    else if(tmp==dist){
                        inserer(liste1,a,b);
                    }
                }
            }
            Bmu *ceBm = liste1->premier;
            if(liste1->nb > 1){
                int _rand = rand()%liste1->nb;
                for(int o=1; o<_rand;o++){
                    ceBm = ceBm->bmusuiv;
                }
            }
            voisinage(ceBm,alpha,ligne,colonne,tabShuffle[j]);
            affecte_val(tabShuffle[j],ceBm);
        }
    }
}


int main(){
    init_vecteur();
    vecteurMoyen();
    vecteurMinimum(0.05);
    vecteurMaximum(0.05);

    init_reseauMap();

    printf("--------------Initialisation-------------\n\n");
    affichage_reseau();

    int ligne =6;
    int colonne =10;
    int nb_iteration = 2000;
    int phase_1= nb_iteration/4;
    int phase_2 = nb_iteration - phase_1;

    printf("\n\n----------Phase 1----------\n\n");

    double alpha = 0.7;
    apprentissage(phase_1,ligne,colonne,alpha);
    affichage_reseau();

    printf("\n\n----------Phase 2----------\n\n");
    alpha = 0.07;
    apprentissage(phase_2,ligne,colonne,alpha);
    affichage_reseau();
    
    printf("\nA correspond à Iris-setosa\n");
    printf("B correspond à Iris-versicolor\n");
    printf("C correspond à Iris-virginica\n");

    free(vecMin);
    free(vecMax);
    free(vecMoyen);

    
    return 0;
}