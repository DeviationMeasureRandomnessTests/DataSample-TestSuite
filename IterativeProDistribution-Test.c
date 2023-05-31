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
/*循环赋值，最后一个参数必须是NULL
 *加入有四个参数a,b,c,NULL，则完成
 * temp = a, a = b, b = c, c= temp的操作
 * 注意，传进的参数值须是一个指针
 *
 * */
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
/*初始化概率数值对链表
 *使链表具有一个表头（不存储数据）
 * */
void InitLinkList(ProbPair **head) {
    *head = (ProbPair *)malloc(sizeof(ProbPair));
    (*head)->h = -1;
    (*head)->p = -1;
    (*head)->next = NULL;
}
/*删除参数Pre所指向的某链表中节点的下一个节点
 *并将后面的节点接上
 * */
void deleteNextOne(ProbPair **Pre) {
    ProbPair *temp = (*Pre)->next;
    (*Pre)->next = temp->next;
    free(temp);
}
/*论文中给出的圈乘运算
 * head:概率分布
 * h:随机变量的取值
 * p:概率值
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
/*向链表中插入取值-概率对
 *
 * */
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
/*两个链表的合并操作，即论文中所说关于概率分布的并运算
 *注意到该操作的结果会使head1链表上的节点全部合并到head0上，即该函数结束后，整个head1链表将不复存在
 *
 * */
void Union(ProbPair *head0, ProbPair * head1) {
    ProbPair *p0, *p1;
    ProbPair *tail;
    p0 = head0;
    p1 = head1;
    while (p1->next != NULL) {
        while (p0->next != NULL && p1->next->h > p0->next->h)
            Next(&p0);
        if (p0->next == NULL) { //追加至其尾
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
/*计算期望
 * */
double Expectation(ProbPair *head) {
    double sum = 0.0;
    ProbPair *temp;
    for (temp = head->next; temp; Next(&temp))
        sum += temp->h * temp->p;
    return sum;
}
/*众数
 * */
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
/*中位数
 * */
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

/*方差
 * */
double Variance(ProbPair *head) {
    double sum = 0.0;
    double e = Expectation(head);
    ProbPair *temp;
    for (temp = head->next; temp; Next(&temp))
        sum +=  temp->p * (temp->h - e) * (temp->h - e);
    return sum;
}

/*计算上分位点
 * 比如sump = 0.1,则计算上10%分位点
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
/*计算下分位点
 * */
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
/*计算head所表示的概率分布在区间[start, end]中的概率
 * */
double Fenwei(ProbPair *head, double start, double end) {
    double sum = 0.0;
    ProbPair *temp;
    for (temp = head->next; temp; Next(&temp)) {
        if (temp->h >= start && temp->h <= end)
            sum += temp->p;
    }
    return sum;
}
/*打印概率分布
 * */
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
/*摧毁链表，并释放其所有节点所占内存*/
void destructLinkList(ProbPair *head) {
    ProbPair *del, *temp;
    for (temp = head; temp->next; ) {
        del = temp->next;
        temp->next = del->next;
        free(del);
    }
    free(head);
}


/*计算长度为M的二元序列一阶离差和的概率分布，在此过程中会计算出长度为M-2,M-4,......的二元序列一阶离差和的概率分布
 *因此返回的是一个概率分布数组, 返回值是G1, G1[M]才表示长度为M的二元序列一阶离差和的概率分布
 * */
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
/*计算长度为M的二元序列二阶离差和的概率分布
 * */
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

/* 计算长度为M的三角形面积和的概率分布
 * 这里参数必须是偶数
 * 注意在M大于0时得到的概率分布A[M]的所有项概率和是0.5，而不是1
 * */
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


/*计算跳跃复杂的概率分布
 * */
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
/*计算长度为M的二元序列的奇序跳跃和与偶序跳跃和的概率分布
 *由于在计算二者的概率分布时，它们会相互用到对方的概率分布，因此索性将二者的概率分布放在同一函数中计算
 *这样一来，会产生两个链表数组，因此通过改变参数pO和pE所指向的链表数组将结果置于其中
 * */
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


/*释放链表数组所占的内存
 *
 * */
void freeDistributions(ProbPair ** G, int M){
	int i; 
    for( i = M & 1; i <= M; i+=2)
        destructLinkList(G[i]); 
    free(G);
}
int main(int argc, char**argv) {
    int i, j, k, M;
   
    /*示例：计算长度为5的二元序列的一阶离差和的概率分布并输出*/
    
    
    M = 436; 
    ProbPair **G1 = dm1ProDistribution(M); 
    PrintLinkList(G1[M]); 
    
    printf("Median: M: %d--%6f\n", M, Median(G1[M]) );

    freeDistributions(G1, M);  
    

    /*示例：计算长度为4的二元序列的二阶离差和的概率分布并输出*/
    
    /*
    M = 4; 
    ProbPair **G2 = dm2ProDistribution(M); 
    PrintLinkList(G2[M]); 
    
    freeDistributions(G2, M); 
    */



    /*示例：计算长度为4的二元序列的奇序跳跃和与偶序跳跃和的概率分布*/
     
    /*
    M = 4; 
    ProbPair **O, **E; 
    ohsEhsProDistribution(M, &O, &E); 
    PrintLinkList(O[M]); 
    PrintLinkList(E[M]); 

    freeDistributions(O, M); 
    freeDistributions(E, M); 
    */

    /*示例：计算长度为4的二元序列的跳跃复杂度的概率分布并输出 */
    
    /*
    M = 4; 
    ProbPair **J = jumpProDistribution(M); 
    PrintLinkList(J[M]); 
    
    freeDistributions(J, M); 
    */


    return 0;
}

