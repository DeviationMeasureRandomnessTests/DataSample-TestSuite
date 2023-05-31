#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<stdarg.h>
#include<math.h>
#include<time.h>
#define endl printf("\n")
#define MALLOC1(A) (A*)malloc(sizeof(A))
struct ProbPair_ {
    double h;
    double p;
    struct ProbPair_ *next;
};
typedef struct ProbPair_ ProbPair;
void Next(ProbPair **p) {
    *p = (*p)->next;
}

void CycleAssign(ProbPair**a, ...) {
    va_list ap;
    va_start(ap, a);
    ProbPair **first, **second;
    ProbPair *temp;
    temp = *a;
    first = a;
    while ((second = va_arg(ap, ProbPair**)) != NULL) {
        *first = *second;
        first = second;
    }
    if (first != a)
        *first = temp;
    va_end(ap);
}
/*Initialize List */
void InitLinkList(ProbPair **head) {
    *head = (ProbPair *)malloc(sizeof(ProbPair));
    (*head)->h = -1;
    (*head)->p = -1;
    (*head)->next = NULL;
}

void deleteNextOne(ProbPair **Pre) {
    ProbPair *temp = (*Pre)->next;
    (*Pre)->next = temp->next;
    free(temp);
}
/*Dot Operation
 * head:Probability Distribution
 * h: Value
 * p: Probability 
 * */
ProbPair * dot(ProbPair *head, double h, double p) {
    ProbPair *c;
    ProbPair *temp, *tempc;
    InitLinkList(&c);
    tempc = c;
    for (temp = head->next; temp != NULL; Next(&temp), Next(&tempc)) {
        ProbPair *newPoint = MALLOC1(ProbPair);
        newPoint->h = temp->h + h;
        newPoint->p = temp->p * p;

        newPoint->next = tempc->next;
        tempc->next = newPoint;

    }
    return c;
}
/*Insert */
void Insert(ProbPair *head, double h, double p) {
    ProbPair * temp;
    for (temp = head; temp->next != NULL && h > temp->next->h; Next(&temp))
        ;
    if (temp->next == NULL) {
        ProbPair * newPoint = (ProbPair*)malloc(sizeof(ProbPair));
        newPoint->h = h;
        newPoint->p = p;
        newPoint->next = temp->next;
        temp->next = newPoint;
    }
    else if (temp->next->h == h) {
        temp->next->p += p;
    }
    else {
        ProbPair * newPoint = (ProbPair*)malloc(sizeof(ProbPair));
        newPoint->h = h;
        newPoint->p = p;
        newPoint->next = temp->next;
        temp->next = newPoint;
    }
}
/*Union */
void Union(ProbPair *head0, ProbPair * head1) {
    ProbPair *p0, *p1;
    ProbPair *tail;
    p0 = head0;
    p1 = head1;
    while (p1->next != NULL) {
        while (p0->next != NULL && p1->next->h > p0->next->h)
            Next(&p0);
        if (p0->next == NULL) { 
            p0->next = p1->next;
            break;
        }
        if (p1->next->h == p0->next->h) {
            p0->next->p += p1->next->p;
            Next(&p0);
            deleteNextOne(&p1);
        }
        else {
            for (tail = p1->next; tail->next != NULL && tail->next->h < p0->next->h; Next(&tail))
                ;
            CycleAssign(&tail->next, &p0->next, &p1->next, NULL);
            
            /*ProbPair *tmp = tail->next;
            tail->next = p0->next;
            p0->next = p1->next;
            p1->next = tmp;*/

            p0 = tail;
        }
    }
    free(head1);
}
/*Calculate Expectation */
double Expectation(ProbPair *head) {
    double sum = 0.0;
    ProbPair *temp;
    for (temp = head->next; temp; Next(&temp))
        sum += temp->h * temp->p;
    return sum;
}
/*Calculate Mode */
double Mode(ProbPair *head) {
    double h_most = head->next->h;
    double p_most = head->next->p;
    ProbPair *temp;
    for (temp = head->next->next; temp; Next(&temp))
        if (temp->p > p_most) {
            h_most = temp->h;
            p_most = temp->p;
        }
    return h_most;
}
/*Calculate Median */
double Median(ProbPair *head) {
    double sum = 0.0;
    ProbPair *temp;
    for (temp = head->next; temp; Next(&temp)) {
        sum += temp->p;
        if (sum >= 0.5)
            break;
    }
    return temp->h;
}

/*Calculate Variance */
double Variance(ProbPair *head) {
    double sum = 0.0;
    double e = Expectation(head);
    ProbPair *temp;
    for (temp = head->next; temp; Next(&temp))
        sum +=  temp->p * (temp->h - e) * (temp->h - e);
    return sum;
}

/*Calculate Upper Percentile, for example, sump = 0.1, 10% Upper Percentile
 * */
double Upper(ProbPair *head, double sump) {
    double sum = 0.0;
    ProbPair *temp;
    for (temp = head->next; temp; Next(&temp)) {
        sum +=  temp->p;
        if (sum >= 1 - sump)
            break;
    }
    return temp->h;
}
/*Calculate Lower Percentile */
double Lower(ProbPair *head, double sump) {
    double sum = 0.0;
    ProbPair *temp;
    for (temp = head->next; temp; Next(&temp)) {
        sum += temp->p;
        if (sum >= sump)
            break;
    }
    return temp->h;
}
/* Calculate  the Probability of head Prob. Distribution within [start, end] */
double Fenwei(ProbPair *head, double start, double end) {
    double sum = 0.0;
    ProbPair *temp;
    for (temp = head->next; temp; Next(&temp)) {
        if (temp->h >= start && temp->h <= end)
            sum += temp->p;
    }
    return sum;
}
/*Print Probability Distribution */
void PrintLinkList(ProbPair *head) {
    int i = 0;
    ProbPair *temp;
    double sum = 0;
    for (temp = head->next; temp; Next(&temp)) {
        printf("(%f, %f) ", temp->h, temp->p);
        sum += temp->p;
        if (++i % 5 == 0 || temp->next == NULL) endl;
    }
    endl;
    printf("total: %d items, the sum of probabilities is %f\n", i, sum);
    printf("----------------------------------------------------------\n");
}
/*Free Memory*/
void destructLinkList(ProbPair *head) {
    ProbPair *del, *temp;
    for (temp = head; temp->next; ) {
        del = temp->next;
        temp->next = del->next;
        free(del);
    }
    free(head);
}


/* Compute the Probability Distribution of Deviation Measure of order 1, given a blocklength M
   Return value G1, G1[M] is the Probability Distribution of Deviation Measure of order 1 with blocklength M */
ProbPair **dm1ProDistribution(int M){
    int i, j, k; 
    ProbPair **G1 = (ProbPair **)malloc(sizeof(ProbPair) * (M+1));
    for(i = M & 1; i <= M; i += 2)
        InitLinkList(&G1[i]);
    
    if(M & 1)//the case that M is odd
        Insert(G1[1], 0.5, 1.0); 
    else     //the case that M is even
        Insert(G1[0], 0, 1.0); 
    int start = (M & 1) + 2; 
    for(i = start; i <= M; i+=2){
        for(j = 1; j <= i/2; j++)
            Union(G1[i], dot(G1[i - 2*j], j*j/2.0, pow(0.5, j))); 
        for(     ; j < i; j++){
            k = 2 * j - i - 1; 
            Insert(G1[i], (2*j*j - k*(k+1)) / 4.0, pow(0.5, j)); 
        }
        k = 2 * j - i - 1; 
        Insert(G1[i], (2*j*j - k*(k+1)) / 4.0, pow(0.5, j-1)); 
    }
    return G1; 
}
/*Compute the Probability Distribution of Deviation Measure of order 2 */
ProbPair ** dm2ProDistribution(int M){
    int i, j, k; 
    ProbPair **G2 = (ProbPair **)malloc(sizeof(ProbPair) * (M+1));
    for(i = M & 1; i <= M; i += 2)
        InitLinkList(&G2[i]);
    
    if(M & 1)//the case that M is odd
        Insert(G2[1], 0.25, 1.0); 
    else     //the case that M is even
        Insert(G2[0], 0, 1.0); 
    int start = (M & 1) + 2; 
    for(i = start; i <= M; i+=2){
        for(j = 1; j <= i/2; j++)
            Union(G2[i], dot(G2[i - 2*j], j*(2*j*j + 1)/12.0, pow(0.5, j))); 
        
        for(     ; j < i; j++){
            k = 2 * j - i - 1; 
            Insert(G2[i], (2*j*(2*j*j+1) - k*(k+1)*(2*k+1)) / 24.0, pow(0.5, j));
        }
        k = 2 * j - i - 1; 
        Insert(G2[i], (2*j*(2*j*j+1) - k*(k+1)*(2*k+1)) / 24.0, pow(0.5, j-1)); 
    }
    return G2; 
}

/* Compute the Probability Distribution of Area, given the blocklength M, where M is even
Note that the sum of probability of all items in A[M] is 0.5, not 1 */
ProbPair ** areaProDistribution(int M){
    if(M & 1){
        printf("M must be even\n"); 
        return NULL; 
    }
    int i, j, k; 
   
    ProbPair **A = (ProbPair **)malloc(sizeof(ProbPair) * (M+1));
    for(i = 0; i <= M; i += 2)
        InitLinkList(&A[i]);
    Insert(A[0], 0, 1.0); 
    for(i = 0; i <= M; i += 2)
        for(j = 1; j <= i/2; j++)
            Union(A[i], dot(A[i - 2*j], j*j/2.0, pow(0.5, j))); 
    return A; 
}


/*Compute the Probability Distribution of Jump complexity */
ProbPair ** jumpProDistribution(int M){
    int i, j, k; 
    ProbPair **J = (ProbPair **)malloc(sizeof(ProbPair) * (M+1));
    for(i = M & 1; i <= M; i += 2)
        InitLinkList(&J[i]);
    
    if(M & 1)//the case that M is odd
        Insert(J[1], 1.0, 1.0); 
    else     //the case that M is even
        Insert(J[0], 0, 1.0); 
    int start = (M & 1) + 2; 
    for(i = start; i <= M; i+=2){
        for(j = 1; j <= i/2; j++)
            Union(J[i], dot(J[i - 2*j], 1.0, pow(0.5, j))); 
        Insert(J[i], 1.0, pow(0.5, i/2) - pow(0.5, i));  
        Insert(J[i], 0, pow(0.5, i)); 
    }
    return J;     
}
/* Compute the Probability Distribution of Odd hop sum, and even hop sum,
where two list arrays(pO & pE) are included*/
void ohsEhsProDistribution(int M, ProbPair ***pO, ProbPair***pE){
    int i, j, k; 
    *pO = (ProbPair **)malloc(sizeof(ProbPair) * (M+1));
    *pE = (ProbPair **)malloc(sizeof(ProbPair) * (M+1));
    ProbPair **O, **E; 
    O = *pO; 
    E = *pE; 
    for(i = M & 1; i <= M; i += 2){
        InitLinkList(&O[i]);
        InitLinkList(&E[i]); 
    }
    if(M & 1){//the case that M is odd
        Insert(O[1], 1.0, 1.0);
        Insert(E[1], 0, 0.0); 
    }
    else{     //the case that M is even
        Insert(O[0], 0, 1.0);
        Insert(E[0], 0, 1.0);
    }
    int start = (M & 1) + 2; 
    for(i = start; i <= M; i += 2){
        for(j = 1; j <= i/2; j++)
            Union(O[i], dot(E[i - 2*j], j, pow(0.5, j))); 
        for(    ; j <= i; j++)
            Insert(O[i], j, pow(0.5, j)); 
        Insert(O[i], 0, pow(0.5, i)); 
    
        for(j = 1; j <= i/2; j++)
            Union(E[i], dot(O[i - 2*j], 0, pow(0.5, j))); 
        Insert(E[i], 0, pow(0.5, i/2)); 
    }

}


/*Free Memory of list array */
void freeDistributions(ProbPair ** G, int M){
	int i; 
    for( i = M & 1; i <= M; i+=2)
        destructLinkList(G[i]); 
    free(G);
}
int main(int argc, char**argv) {
    int i, j, k, M;
   
    
    
    
    M = 436; 
    ProbPair **G1 = dm1ProDistribution(M); 
    PrintLinkList(G1[M]); 
    
    printf("Median: M: %d--%6f\n", M, Median(G1[M]) );

    freeDistributions(G1, M);  
    

    /*Example£º Calculate the Probability Distribution of Deviation Measure of order 2, for M=4*/    
    /*
    M = 4; 
    ProbPair **G2 = dm2ProDistribution(M); 
    PrintLinkList(G2[M]); 
    
    freeDistributions(G2, M); 
    */


     /*Example£º Calculate the Probability Distribution of Odd hop sum and Even hop sum, for M=4*/        
    /*
    M = 4; 
    ProbPair **O, **E; 
    ohsEhsProDistribution(M, &O, &E); 
    PrintLinkList(O[M]); 
    PrintLinkList(E[M]); 

    freeDistributions(O, M); 
    freeDistributions(E, M); 
    */

     /*Example£º Calculate the Probability Distribution of jump complexity, for M=4*/
        
    /*
    M = 4; 
    ProbPair **J = jumpProDistribution(M); 
    PrintLinkList(J[M]); 
    
    freeDistributions(J, M); 
    */


    return 0;
}

