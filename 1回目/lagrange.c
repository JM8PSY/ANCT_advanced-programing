#include <stdio.h>
#include <stdlib.h>

#define STEPS 10000

struct point {
    double x,y;
};

double lagrange(int,int,double,double []);

int main(){
    int i,j;
    int n;
    double t;
    double *tt;

    double lag;

    FILE *input, *output;

    struct point *points;
    
    struct point p;

    input = fopen("points.txt","r");
    output = fopen("output.txt","w");

    printf("次数を入力してください: ");
    scanf("%d",&n);

    points = (struct point *)malloc((n+1)*sizeof(struct point));
    tt = (double *)malloc((n+1)*sizeof(double));

    for(i=0;i<=n;i++){
        fscanf(input,"%lf %lf",&points[i].x,&points[i].y);

        tt[i] = points[i].x;
    }

    // 修正後のループ部分
    for(i=0; i<=STEPS; i++){
        // 未初期化のjではなく，外側のループ変数iを使用する
        t = i * (tt[n] - tt[0]) / STEPS + tt[0];
        p.y = 0.0;
        
        // 内側のループ変数をjに変更（変数iの重複を回避）
        for(j=0; j<=n; j++){
            lag = lagrange(j, n, t, tt);
            p.y += lag * points[j].y;
        }
        // x座標としてt，y座標としてp.yを出力する
        fprintf(output, "%lf %lf\n", t, p.y);
    }

    free(points);
    free(tt);

    fclose(input);
    fclose(output);

    return 0;
}

double lagrange(int i,int n,double t,double tt[]){
    int j;
    double lag = 1.0;

    for(j=0;j<=n;j++){
        if(i!=j){
            lag *= (t - tt[j]) / (tt[i] - tt[j]);
        }
    }

    return lag;
}