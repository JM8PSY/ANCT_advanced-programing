#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NOS 1000

struct point{
    double x,y;
};

void compute_v(struct point *,double *,int);
void elimination(double **,double *,double *,int);
void substitution(double **,double *,double *,int);

int main(void){
    int i,j,n;
    double a,b,c,d,h,dy,t;
    double *v;
    struct point *points;
    struct point p;

    FILE *input, *output;

    printf("次数を入力してください: ");
    if(scanf("%d",&n)!=1 || n<1){
        fprintf(stderr,"次数は1以上の整数で入力してください\n");
        return 1;
    }

    input = fopen("points.txt","r");
    if(input==NULL){
        fprintf(stderr,"points.txt を開けません\n");
        return 1;
    }
    output = fopen("curve.txt","w");
    if(output==NULL){
        fprintf(stderr,"curve.txt を開けません\n");
        fclose(input);
        return 1;
    }

    points = (struct point *)malloc((n+1)*sizeof(struct point));
    v      = (double *)malloc((n+1)*sizeof(double));
    if(points==NULL || v==NULL){
        fprintf(stderr,"メモリ確保に失敗しました\n");
        return 1;
    }

    for(i=0;i<=n;i++){
        if(fscanf(input,"%lf %lf",&points[i].x,&points[i].y)!=2){
            fprintf(stderr,"点データが足りません(%d 点必要)\n",n+1);
            return 1;
        }
    }

    /* 自然スプライン:両端の2次係数は0 */
    v[0] = 0.0;
    v[n] = 0.0;

    if(n>=2){
        compute_v(points,&v[1],n-1);
    }

    for(i=0;i<n;i++){
        h  = points[i+1].x - points[i].x;
        dy = points[i+1].y - points[i].y;

        a = points[i].y;                          /* 定数項 */
        b = v[i];                                 /* 2次の係数 */
        c = dy/h - h*(2.0*v[i]+v[i+1])/3.0;       /* 1次の係数 */
        d = (v[i+1]-v[i])/(3.0*h);                /* 3次の係数 */

        for(j=0;j<NOS;j++){
            t = j/(double)NOS*h;
            p.x = points[i].x + t;
            p.y = a + c*t + b*t*t + d*t*t*t;
            fprintf(output,"%lf %lf\n",p.x,p.y);
        }
    }
    /* 最後の点 */
    fprintf(output,"%lf %lf\n",points[n].x,points[n].y);

    free(points);
    free(v);
    fclose(input);
    fclose(output);

    return 0;
}

/* 内部節点 1..n の2次係数 v を求める(n は未知数の個数) */
void compute_v(struct point *points,double *v,int n){
    int i,k;
    double **A,*b;
    double h0,h1;

    A = (double **)malloc(n*sizeof(double *));
    for(i=0;i<n;i++){
        A[i] = (double *)calloc(n,sizeof(double));   /* ゼロ初期化 */
    }
    b = (double *)malloc(n*sizeof(double));

    for(i=0;i<n;i++){
        k  = i+1;                                  /* 対応する節点番号 */
        h0 = points[k].x   - points[k-1].x;
        h1 = points[k+1].x - points[k].x;

        A[i][i] = 2.0*(h0+h1);
        if(i>0)   A[i][i-1] = h0;
        if(i<n-1) A[i][i+1] = h1;

        b[i] = 3.0*((points[k+1].y-points[k].y)/h1
                  - (points[k].y-points[k-1].y)/h0);
    }

    elimination(A,b,v,n);

    for(i=0;i<n;i++){
        free(A[i]);
    }
    free(A);
    free(b);
}

void elimination(double **A,double *b,double *v,int n){
    int i,j,k;
    double m;

    for(k=0;k<n-1;k++){
        for(i=k+1;i<n;i++){
            m = A[i][k]/A[k][k];
            for(j=k;j<n;j++){
                A[i][j] -= m*A[k][j];
            }
            b[i] -= m*b[k];
        }
    }

    substitution(A,b,v,n);
}

void substitution(double **A,double *b,double *v,int n){
    int i,j;

    for(i=n-1;i>=0;i--){
        v[i] = b[i];
        for(j=i+1;j<n;j++){
            v[i] -= A[i][j]*v[j];
        }
        v[i] /= A[i][i];
    }
}