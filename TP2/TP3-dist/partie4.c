#define N 12

int main(){
double A[N], B[N];

for(int i=0 ; i<N; i++){
    if(i&0x1) { // i impair
        B[i]=B[i]+A[i] ;
    } else { // i pair
        B[i]=B[i]-A[i] ;
    }
}

return 0;
}