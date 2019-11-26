## **StarHai Game Theory总结**

## **Nim Game**

### **Nim**

**问题**：共有N堆石子，编号1..n，第i堆中有个a[i]个石子。

每一次操作Alice和Bob可以从任意一堆石子中取出任意数量的石子，至少取一颗，至多取出这一堆剩下的所有石子。

**结论**：对于一个局面，当且仅当a[1] xor a[2] xor ...xor a[n]=0时，该局面为P局面，即必败局面。

**证明**：二进制位证明即可。

### **Moore’s Nim**

**问题**：n堆石子，每次从不超过k堆中取任意多个石子，最后不能取的人失败。

**结论**：这是一个nim游戏的变形：把n堆石子的石子数用二进制表示，统计每个二进制位上1的个数，若每一位上1的个数mod(k+1)全部为0，则必败，否则必胜。(先手)

**证明**：分类讨论N/P状态。

### **Staircase N****im**

**问题**：在阶梯上进行，每层有若干个石子，每次可以选择任意层的任意个石子将其移动到该层的下一层。最后不能操作的人输。

**结论**：在奇数堆的石子做Nim。

**证明**：阶梯博弈经过转换可以变为Nim.把所有奇数阶梯看成N堆石子做nim。把石子从奇数堆移动到偶数堆可以理解为拿走石子，就相当于几个奇数堆的石子在做Nim。

### **New** **Nim**

**问题**：在第一个回合中，第一个游戏者可以直接拿走若干个整堆的火柴。可以一堆都不拿，但不可以全部拿走。第二回合也一样，第二个游戏者也有这样一次机会。从第三个回合（又轮到第一个游戏者）开始，规则和Nim游戏一样。如果你先拿，怎样才能保证获胜？如果可以获胜的话，还要让第一回合拿的火柴总数尽量小。

**结论**：为使后手必败，先手留给后手的必然是若干线性无关的数字，否则后手可以留下一个异或和为零的非空子集使得先手必败，故问题转化为拿走和最小的数字使得留下的数线性无关，即留下和最大的线性基，这样拿走的数量显然最少，找到和最大的线性基只需贪心的把数字从大到小加入到基中即可.

**证明**:证明需用到拟阵,但结论的正确性是比较显然的.

### **A****nti-****N****im**

**问题**：正常的nim游戏是取走最后一颗的人获胜，而反nim游戏是取走最后一颗的人输。

**结论**：一个状态为必胜态，当且仅当：

　　1）所有堆的石子个数为1，且NIM_sum(xor和)=0

　　2）至少有一堆的石子个数大于1，且 NIM_sum(xor和) ≠ 0

![img](file:///C:\Users\Lenovo\AppData\Local\Temp\ksohtml39104\wps1.jpg) 

### **Lasker’s Nim**

**问题**：Alice和Bob轮流取石子，每一次可以从任意一堆中拿走任意个石子，也可以将一堆石子分为两个小堆。先拿完者获胜。

结论：if(x%4==0) sg[x]=x-1;if(x%4==1||x%4==2) sg[x]=x;if(x%4==3) sg[x] = x+1。**这种问题一般要对sg值打表。(对于****先拿完者获胜**：**sg值为零先手必败，否则后手必胜)**

## **Wythoff's game**

**问题**：两堆石子，每次可以取一堆或两堆，从两堆中取得时候个数必须相同，先取完的获胜。

**结论**： ak =[k *（1+√5）/2]，bk= ak + k  （k=0，1，2，…,n 方括号表示取整函数)

**证明**:

这种情况下是颇为复杂的。我们用（ak，bk）（ak ≤ bk ,k=0，1，2，…,n)表示两堆物品的数量并称其为局势，如果甲面对（0，0），那么甲已经输了，这种局势我们称为奇异局势。前几个奇异局势是：（0，0）、（1，2）、（3，5）、（4，7）、（6，10）、（8，13）、（9，15）、（11，18）、（12，20）。

​    可以看出,a0=b0=0,ak是未在前面出现过的最小自然数,而 bk= ak + k，奇异局势有

如下三条性质：

​    1.任何自然数都包含在一个且仅有一个奇异局势中。

​    由于ak是未在前面出现过的最小自然数，所以有ak > ak-1 ，而 bk= ak + k > ak-1 + k-1 = bk-1 > ak-1 。所以性质1。成立。

​    2.任意操作都可将奇异局势变为非奇异局势。

​    事实上，若只改变奇异局势（ak，bk）的某一个分量，那么另一个分量不可能在其他奇异局势中，所以必然是非奇异局势。如果使（ak，bk）的两个分量同时减少，则由于其差不变，且不可能是其他奇异局势的差，因此也是非奇异局势。

​    3.采用适当的方法，可以将非奇异局势变为奇异局势。  

​    从如上性质可知，两个人如果都采用正确操作，那么面对非奇异局势，先拿者必胜；反之，则后拿者取胜。

4.（Betty 定理）：如果存在正无理数 A, B 满足 1/A + 1/B = 1，那么集合 P = { [At], t ∈ Z+}、Q = { [Bt], t ∈ Z+} 恰为集合 Z+ 的一个划分，即：P ∪ Q = Z+，P ∩ Q = ø。

5.上述矩阵中每一行第一列的数为 [Φi]，第二列的数为 [(Φ + 1)i]，其中 Φ = (sqrt(5) + 1) / 2 为黄金分割比。

## **Bash Game**

问题:只有一堆石子共n个。每次从最少取1个，最多取m个，最后取光的人取胜。

结论：如果n=(m+1)*k+s (s!=0) 那么先手一定必胜。

证明：因为第一次取走s个，接下来无论对手怎么取，我们都能保证取到所有(m+1)倍数的点，那么循环下去一定能取到最后一个。

 

## **Sprague-Grundy Theorem**

**Sprague-Grundy Theorem：**

**SG(x)=mex{SG(y) | y是x的后继}{SG(x):表示当前拿多少为必胜}**

g(G)=g(G1)^g(G2)^...^g(Gn)。也就是说，游戏的和的SG函数值是它的所有子游戏的SG函数值的异或。

 

### **Fibonacci’s Game**

**任何正整数可以表示为若干个不连续的Fibonacci数之和。**

问题：有一堆个数为n的石子，游戏双方轮流取石子，满足：

1)先手不能在第一次把所有的石子取完；

2)之后每次可以取的石子数介于1到对手刚取的石子数的2倍之间（包含1和对手刚取的石子数的2倍）。

结论：先手胜当且仅当n不是Fibonacci数。换句话说，必败态构成Fibonacci数列。

### **Game On The Tree**

#### **树链博弈：**

给定一棵 n 个点的树，其中 1 号结点是根，每个结点要么是黑色要么是白色
现在小 Bo 和小 Biao 要进行博弈，他们两轮流操作，每次选择一个黑色的结点将它变白，之后可以选择任意多个(可以不选)该点的祖先(不包含自己)，然后将这些点的颜色翻转，不能进行操作的人输
由于小 Bo 猜拳经常输给小 Biao，他想在这个游戏上扳回一城，现在他想问你给定了一个初始局面，是先手必胜还是后手必胜。

题解：每层的黑点数为偶数的时候，为先手必败态。首先，没有黑点，都是0，肯定是先手必败态。

 

 

####  **Bamboo Stalks**

 	作为GH游戏的介绍，我们先研究下面的Figure 6.1。n条线段的bamboo stalks游戏是具有n条边的线形图。一步合法的操作是移除任意一条边，玩家轮流进行操作，最后一个进行操作的玩家获胜。n条线段的bamboo stalks游戏能够移动到任意更小线段数(0到n-1)的bamboo stalks游戏局面当中。所以n条线段的bamboo stalks游戏是等同于nim游戏中其中拥有n个石子的一堆。玩一组bamboo stalks游戏就相当于玩nim游戏。

![img](file:///C:\Users\Lenovo\AppData\Local\Temp\ksohtml39104\wps2.png) 

例如，左边的三根竹竿构成的“森林”相当于具有石子数分别为3、4、5三堆石子的nim游戏。就我们所知，3^4^5=2，这是一个能够移动到P局面的N局面，办法是通过取走三根线段的竹竿上的第二根线段，留下一根。而结果变成右边的竹竿分布，而此时的SG值是0，是P局面。

#### **Green Hackenbush on Trees**

**Colon Principle****：当树枝在一个顶点上时，用一个非树枝的杆的长度来替代，相当于他们的****n****异或之和。**

通过bamboo stalks游戏，我们知道GH游戏是换了个形式的nim游戏而已。可是如果我们允许比bamboo stalks游戏更多结构呢？让我们看下在Figure 6.2中由三棵根树组成的森林。根树是一种图，带有一个最高的节点，叫做根，这个根到任意一个其他节点的路径都是独一无二的。实质上就是说这个图不含有圈。

![img](file:///C:\Users\Lenovo\AppData\Local\Temp\ksohtml39104\wps3.png) 

一次合法操作是移除任意一条不与地面相连的线段，此时次线段以上的子树都会被移除。既然这个游戏是公平的，而且根据我们学过的nim游戏，这样的树相当于nim游戏的一些堆，或者说是bamboo stalks游戏（到了这里明白bamboo stalks游戏与nim游戏的单堆是等价的）。这个问题就是寻找每棵树的SG值。

​          **这里我们要用上一个原理，叫做****Colon Principle****：当树枝在一个顶点上时，用一个非树枝的杆的长度来替代，相当于他们的****n****异或之和。**

​         让我们看看这条原理是如何在Figure 6.2中寻找与左树等价的竹竿。这里有两个节点拥有树枝。较高的那个节点拥有两个树枝，每个树枝上有两个节点。1^1=0，所以这两个树枝可以用一个带有0个节点树枝来代替，也就是说可以把这两个树枝给去掉。那就剩下了一个Y形树了，因为同样道理这个Y形树的两个树枝也要被去掉。此时就剩下了一个线段数为1的竹竿游戏模型了。

​          看Figure 6.3，是Figure 6.2中第二棵树的处理办法，最后可以得到线段数为8的竹竿游戏。第三棵树也可以同样处理，结果是线段数为4的竹竿游戏。（注意，这里所指的竹竿游戏都实质上是nim游戏中的单堆石子）

![img](file:///C:\Users\Lenovo\AppData\Local\Temp\ksohtml39104\wps4.png) 

现在我们可以计算一下图6.2三棵树的sg值，也就是1^8^4=13.既然这个SG值不为0，那么就是一个N局面，先手必胜。问题是要怎么找到胜利方法。很明显这里有一个必胜移动，通过将第二棵树进行操作使得它的SG值为5.可是我们要找哪条边呢？

​          Figure 6.3的最后一个树长度是8，因为它的前一棵树的三个树枝长度分别是3,2,6，异或值为3^2^6=7，用长度为7的竹竿代替三个树枝后，树的长度就是根加上竹竿长度，即1+7=8。为了最后SG达到5，即使树的长度为5，我们要用长度为4的竹竿来替代那三个树枝。因为2^6=4，所以我们只要将最左边的树枝去掉就行了，当然我们也可以将树枝改动成3^1^6=4.

修剪树的方法用冒号给出了，把所有的树化简为一个单一的竹竿。一个从最高的树枝开始，然后用原理归纳往下到根部。我们现在**展示这个原理对于含圈和多重根边的图同样适用。**

 

##### **变形一:边权大于1**

LightOJ1355 题目大意：
给一棵带边权的树，两个人分别给边涂色。边权代表了这条边可以涂色的次数。如果一条边的涂色次数没有用完，那么可以涂他子树的边。无法涂色者为负。

题解：green博弈变形，对于都是1的就是green博弈SG[u]^=SG[v]+1;

对于大于1的边，偶数对其没有贡献，奇数有贡献，SG[u]^= SG[v]^(val[v]%2);

 

 

### **Green Hackenbush on** **general rooted graphs**

  **The Fusion Principle****：任何环内的节点可以融合成一点而不会改变图的****sg****值。（下面我们称它为融合原则）**

​        **融合原则允许我们把任意一个根图简化为一个等效的可以通过冒号原则（即****Colon Principle****）简化为竹竿的树。**

 

![img](file:///C:\Users\Lenovo\AppData\Local\Temp\ksohtml39104\wps5.png) 

同样，上面的三个图，每个图都相当于nim游戏的一个堆，三个图组成了一个nim游戏。接下来我们要找到这些图等价的nim堆，方便我们解决问题。这要用到融合原则。我们可以把两个相邻的节点合成一个节点，并把它们之间的边弯曲，变成一个圈。一个圈是把自己作为边的另一端的一种边。比如Figure 6.5的最右边那个杂戏表演者的头就是一个圈。在GH游戏中，一个圈是可以被一个叶子（一条没有任何树枝与它相连的边）所代替的，见Figure 6.6中第三幅图到第四幅图的转化。

  **The Fusion Principle****：任何环内的节点可以融合成一点而不会改变图的****sg****值。（下面我们称它为融合原则）**

​        **融合原则允许我们把任意一个根图简化为一个等效的可以通过冒号原则（即****Colon Principle****）简化为竹竿的树。**

如下图左部分所示的一个门，在地板上的两个节点是同样的节点来的（记住地板相当于一个单独的节点），所以实际上是一个有一个节点与地板相连的三角形，即第二幅图。融合原则告诉我们，这相当于一个单独的节点有三个圈与它相连。所以造就了第三幅到第四幅的转变，过程是把任意两点收缩成一个圈，3个点两两收缩便可得到三个圈。每个圈又相当于长度为1的竹竿，它们的异或和还是长度为1的竹竿。

![img](file:///C:\Users\Lenovo\AppData\Local\Temp\ksohtml39104\wps6.png) 

 我们会发现，拥有奇数条边的环可简化为一条边，偶数条边的环可简化为一个节点。例如，在Figure 6.5中的第二幅图圣诞树中的有四条边的环，会缩减成一个节点，所以这圣诞树最后会简化为一个长度为1的竹竿。相似的，房子上的烟囱变成一个单独的节点，右边的窗户变成一个点，继续下去，就可以看出房子的SG值为3。

![img](file:///C:\Users\Lenovo\AppData\Local\Temp\ksohtml39104\wps7.png) 

 

 

 

 

 

### **SG函数的一些题目**

LightOJ1199

题意：有n堆石子(1<=n<=100)，每一堆分别有xi个石子(1<=xi<=10000)， 一次操作可以使一堆石子变成两堆数目不相等的石子， 最后不能操作的算输，问先手胜还是后手胜。

对于每一个数下一个拆分为 (1,n-1),(2,n-2) ...  个.SG函数推到即可；

```c++
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
const int maxn=10010;
int T,n,SG[maxn],x,ans,vis[maxn]; 
void getSG(int n)
{
    for(int i=1;i<=n;++i)
    {
        memset(vis,0,sizeof vis);
        for(int j=1;j*2<i;++j) if(i!=j*2) vis[SG[j]^SG[i-j]]=1;
        for(int j=0;j<maxn;++j){ if(!vis[j]) { SG[i]=j;break; } }
    }
}
int main()
{
    scanf("%d",&T);
    getSG(maxn);
    for(int cas=1;cas<=T;++cas)
    {
        scanf("%d",&n);ans=0;
        for(int i=1;i<=n;++i)
        {
            scanf("%d",&x);
            ans^=SG[x];
        }
        if(ans) printf("Case %d: Alice\n",cas);
        else printf("Case %d: Bob\n",cas);
    }
    return 0;    
}
```

 

LightOJ1229

题目意思就是给你一行字符串由 '.'和'X'组成，然后两个人交替将一个‘.’

变成'X'.如果某个人先形成连续的3个‘X’，则这个人就取得胜利。问先手必胜的位置是否存在，

如果存在，有多少个，并依次输出其位置；这题肯定是枚举每一个位置，判断是否可以胜利，

对于SG[x]表示长度为x的‘.’区间的SG值，然后对于每一个位置（不是'X'的位置），判断其由

‘.’变成'X'之后是否可以形成连续3个'X'的必胜状态，是否会形成.XX或XX.或X.X的必败状态；如果都不是再去枚举每一个区间的SG值,再将其异或ans,就得到这一个位置的SG值，为0必胜，不为零必败；

 

```c++
#include<iostream>
#include<cstring>
#include<cstdio>
#include<vector>
using namespace std;
typedef long long ll;
#define clr(a,val) memset(a,val,sizeof(a))
const int maxn=210;
int T,SG[maxn],len;
vector<int> v;
char s[maxn],s1[maxn];
int getSG(int m)
{
    if(m<0) return 0;
    if(SG[m]!=-1) return SG[m];
    bool vis[maxn];clr(vis,0);
    for(int i=1;i<=m;++i) vis[getSG(i-3)^getSG(m-i-2)]=1;
    int t=0;
    while(vis[t]) ++t;
    return SG[m]=t; 
}

bool check(int x)
{
    strcpy(s1,s);
    if(s1[x]=='X') return 0;
    s1[x]='X';
    for(int i=0;i<len-2;++i) {if(s1[i]=='X'&&s1[i+1]=='X'&&s1[i+2]=='X') return 1;}
    for(int i=0;i<len-1;++i) {if(s1[i]=='X'&&s1[i+1]=='X') return 0;}
    for(int i=0;i<len-2;++i) {if(s1[i]=='X'&&s1[i+2]=='X') return 0;}
    int j=-1,f=0,ans=0;
    for(int i=0;i<len;++i)
    {
        if(s1[i]=='X')
        {
            if(f) ans^=getSG(i-j-5);
            else ans^=getSG(i-j-3),f=1;    
            j=i;
        }
    }
    ans^=getSG(len-j-3);
    return ans==0;
}

int main()
{
    scanf("%d",&T);
    memset(SG,-1,sizeof SG);
    for(int cas=1;cas<=T;++cas)
    {
        scanf("%s",s);
        v.clear(); 
        len=strlen(s);
        for(int i=0;i<len;++i)
        {
            if(check(i)) v.push_back(i+1);
        }    
        printf("Case %d:",cas);
        if(v.size())
        {
            for(int i=0;i<v.size();++i) printf(" %d",v[i]);
            puts("");
        }
        else printf(" 0\n");
    }    
    
    return 0;
}
```

 

LightOJ1344

现在有几串手镯，每个手镯上都有一些珠子，每个珠子都有权值，Aladdin 和 Genie 轮流取珠子，当取了一颗珠子后，这个手镯上所有权值大于等于这颗珠子权值的珠子，都要被删去。因此，这个手镯就会变成新的几个手镯（因为被切割了）。

如，5-1-7-2-4-5-3 这串手镯，选第一颗珠子，权值为5，因此，5，7，5 就要被删去。手镯变成了新的 1，2-4，3 三个手镯。 最后谁不能拿，谁输；

题解：每个手镯都是独立的，因此可以异或每个手镯的 sg 值求解。而每个手镯，可以在某些结点被取走珠子，变成新的几段，任意一种情况的 sg 值，便是新分成的几段的 sg 值异或起来。再将每一种情况的 sg 值，记录在 vis 中，查找没出现过的最小的值，便是这个手镯的 sg 值。

 

```C++
#include<bits/stdc++.h>
using namespace std;
const int maxn = 55;
int sg[maxn][maxn]; 
int arr[maxn][maxn], num[maxn];
int sgtmp[maxn]; 
int ret[maxn][maxn]; 

struct node {
    int x, y;
} output[maxn * maxn]; 
bool operator==(node a,node b){return a.x==b.x&&arr[a.x][a.y]==arr[b.x][b.y];}
bool cmp(node a,node b)
{
    if(a.x==b.x) return arr[a.x][a.y]<arr[b.x][b.y];
    return a.x<b.x;
}

int dfs(int now,int l,int r)
 {
    if(l>r) return 0;
    if(l==r) return sg[l][r]=1; 
    if(sg[l][r]!=-1) return sg[l][r]; 

    int vis[maxn]; 
    memset(vis,0,sizeof(vis));

    for(int i=l;i<=r;++i) 
    {
        int tmp=0,last=-1;

        for(int j=l;j<=r;++j) 
        {
            if(arr[now][j]<arr[now][i]) 
            {
                last = j;
                break;
            }
        }
        for(int j = last + 1; j <= r && last != -1; ++j) 
        {
            if(arr[now][j] >= arr[now][i]) 
            {
                tmp ^= dfs(now, last, j - 1);
                last = -1;
                for (int k = j + 1; k <= r; ++k) 
                {
                    if (arr[now][k] < arr[now][i]) 
                    {
                        last = j = k;
                        break;
                    }
                }
            }
        }
        if(last != -1) tmp ^= dfs(now, last, r);
        vis[tmp] = 1;
        if (l == 1 && r == num[now]) ret[now][i] = tmp;
    }
    for(int i = 0;;++i) 
    {
        if(vis[i]==0) 
        {
            sg[l][r]=i;
            return i;
        }
    }
}
int main() 
{
    int t, k, ca = 1;
    scanf("%d", &t);
    while (t--) 
    {
        int ans = 0;
        scanf("%d", &k);
        for(int i = 1; i <= k; ++i) 
        {
            memset(sg, -1, sizeof(sg));
            scanf("%d", &num[i]);
            for(int j = 1; j <= num[i]; ++j) scanf("%d", &arr[i][j]);
            sgtmp[i] = dfs(i, 1, num[i]);
            ans ^= sgtmp[i];
        }
        if (ans == 0) printf("Case %d: Genie\n", ca++);
        else 
{
            int cnt = 0;
            printf("Case %d: Aladdin\n", ca++);
            memset(output, 0, sizeof(output));
            for (int i = 1; i <= k; ++i)
                for (int j=1;j<=num[i];++j){if((ans^sgtmp[i]^ret[i][j])==0) output[cnt++] = {i, j};}
            sort(output,output+cnt,cmp);
            cnt=(int)(unique(output,output+cnt)-output);
            for(int i = 0; i < cnt; ++i) 
                printf("(%d %d)", output[i].x, arr[output[i].x][output[i].y]);
            printf("\n");
        }
    }
    return 0;
}
```

 