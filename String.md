## **StarHai字符串模板**

 

[TOC]



### **KMP**

```c++
void GetFail(char *P,int *Next)//得到失配函数
{
	int len = strlen(P);
	Next[0]= 0, Next[1]=0;	
	for(int i=1;i<len;i++)
	{
		int j=Next[i];		
		while(j && P[i]!=P[j]) j=Next[j];		
		Next[i+1] = P[i]==P[j]? j+1:0; 
	}	
}

void Find(char *T,char *P,int *Next)
{//查找字符串在文本串出现的位置
	int len1 = strlen(T), len2 = strlen(P);	
	GetFail(P,Next);	
	int j = 0;	
	for(int i=0;i<len1;i++)	
{
		while(j && P[j]!=T[i])  j=Next[j];		
		if(P[j]==T[i])  j++;		
		if(j==len2)  cout<<i-len2+1<<endl;						
	}	
}
```



相似是指长度相同，且如果短字符串中两个位置的字符相同则子串中也相同，不同则也不同，比如：qtjj和well。现在Qt有一个短字符串S，一个长字符串A，他想知道在这个长字符串中有几个子串与短字符串相似。

```C++
#include<bits/stdc++.h>
using namespace std;
const int maxn=2e6+10;
int n;
char s[maxn],p[maxn];
int last[505],tot,la,lb;
int a[maxn],b[maxn],fail[maxn];
bool ans[maxn];

inline bool check(int x,int y)
{
	if(b[x]) return b[x]==a[y];
	else return (a[y]==0)||(a[y]>=x);
}

inline bool check2(int x,int y)
{
	if(b[x]) return b[x]==b[y];
	else return (b[y]==0)||(b[y]>=x);
}

inline void work()//类似KMP 
{
	for(int i=2;i<=lb;++i)
	{
		int j=fail[i-1];
		while(j&&!check2(j+1,i)) j=fail[j];
		if(check2(j+1,i)) fail[i]=j+1;
	}
	int j=0;
	for(int i=1;i<=la;++i)
	{
		while(j&&!check(j+1,i)) j=fail[j];
		if(check(j+1,i)) ++j;
		if(j==lb)
		{
			++tot;
			ans[i-lb+1]=true;
			j=fail[j];
		}
	}
}
 
int main()
{
	scanf("%s",p+1);
	scanf("%s",s+1);
	la=strlen(s+1),lb=strlen(p+1);
	for(int i=0;i<=300;++i) last[i]=0;
	for(int i=1;i<=la;++i)
	{
		if(last[s[i]]) a[i]=i-last[s[i]];
		last[s[i]]=i;
	}
	for(int i=0;i<=300;++i) last[i]=0;
	for(int i=1;i<=lb;++i)
	{
		if(last[p[i]]) b[i]=i-last[p[i]];
		last[p[i]]=i;
	}
	work();
	printf("%d\n",tot);
    return 0;
}
/*
aba
abababa
*/ 
```

 

### **EXKMP**

(求s2与s1的每一个后缀的最大公共前缀)

const int maxn=100010; //字符串长度最大值  

int next[maxn],ex[maxn]; //ex数组即为extend数组  

//预处理计算next数组 

辅助数组next[i]表示s2[i,m-1]和s2的最长公共前缀长度， 

```C++
void GETNEXT(char *str)  
{  
    int i=0,j,po,len=strlen(str);  
    next[0]=len;
    while(str[i]==str[i+1]&&i+1<len) i++;  
	next[1]=i;  
	po=1;    
	for(i=2;i<len;i++)  
    {  
        if(next[i-po]+i<next[po]+po)
        next[i]=next[i-po];  
        else {  
            j=next[po]+po-i;  
            if(j<0) j=0;  
            while(i+j<len&&str[j]==str[j+i])   j++;  
            next[i]=j;  po=i;
        }  
    }  
}  

void EXKMP(char *s1,char *s2)  //计算extend数组  
{  
    int i=0,j,po,len=strlen(s1),l2=strlen(s2);  
    GETNEXT(s2);//计算子串的next数组  
    while(s1[i]==s2[i]&&i<l2&&i<len)//计算ex[0]  
    i++;  ex[0]=i;  
    po=0;//初始化po的位置  
    for(i=1;i<len;i++)  
    {  
        if(next[i-po]+i<ex[po]+po)//第一种情况，直接可以得到ex[i]的值  
        ex[i]=next[i-po];  
        else//第二种情况，要继续匹配才能得到ex[i]的值  
        {  
            j=ex[po]+po-i;  
            if(j<0)j=0;//如果i>ex[po]+po则要从头开始匹配  
            while(i+j<len&&j<l2&&s1[j+i]==s2[j])   j++;  
            ex[i]=j;  
            po=i;//更新po的位置  
        }  
    }  
} 

```

一个数组s[n]={s1,s2..sn},定义h[x]:表示在s[n]后面添加一个数字x后，新增的子串的数量.求：

⨁(c=1,m)(h(c)⋅3cmod(109+7)).

 

```C++
#include<bits/stdc++.h>
using namespace std;
#define mod 1000000007
const int N=1e6+10;
int x[N],y[N],w[N],n,m,nxt[N],ex[N];
void pre_exkmp()
{
    nxt[0]=n;
    int j=0;
    while(j+1<n&&x[j]==x[j+1])j++;
    nxt[1]=j;
    int k=1;
    for(int i=2;i<n;i++)
	{
        int p=nxt[k]+k,t=nxt[i-k];
        if(i+t<p)nxt[i]=t;
        else 
		{
            j=p-i<0?0:p-i;
            while(i+j<n&&x[i+j]==x[j])j++;
            nxt[i]=j;k=i;
        }
    }
}
int main()
{
    while(~scanf("%d%d",&n,&m))
	{
        for(int i=0;i<=n-1;++i) scanf("%d",&x[n-i-1]);
        pre_exkmp();
        for(int i=1;i<=m;++i) w[i]=-1;
        nxt[n]=0;
        for(int i=0;i<=n-1;++i) w[x[i]]=max(nxt[i+1],w[x[i]]);
        int ans=0,p=1;
        for(int i=1;i<=m;++i)
		{
            p=1ll*p*3%mod;
            ans^=1ll*p*(n-w[i])%mod;
        }
        printf("%d\n",ans);
    }
    return 0;
}
```



### **Trie**

(HDU1251查找以某个串为前缀的串的数量)

```C++
struct Trie{
	int ch[maxnode][sidma_size];
	int pos;
	int val[maxnode];
	int num[maxnode];
	trie()
	{
		memset(ch[0],0,sizeof(ch[0]));
		memset(num,0,sizeof(num));
		pos=1;
	} 
	int idx(char chr) { return chr-'a'; }
	void Insert(char *s,int v)
	{
			int u=0,n=strlen(s);
			for(int i=0;i<n;i++)
			{
				int c=idx(s[i]);
				if(!ch[u][c])
				{
					memset(ch[pos],0,sizeof(ch[pos]));
					val[pos]=0;
					ch[u][c]=pos++;
				}
				u=ch[u][c];
				num[u]++;
			}
			val[pos]=v;
	}
	int Query(char s)//	
	{
		int n=strlen(s),u=0;
		for(int i=0;i<n;i++)
		{
			int c=idx(s[i]);
			if(!ch[u][c]) break;
			u=ch[u][c];
		}
		return num[u];
	}
} Tree;

```

n个集合，你要进行m个操作。总共有3种操作。第一种，合并两个集合x和y。第二张，把特定的集合里面所有的数字加一。第三种，询问在某个集合里面，对于所有数字对2的k次方取模后，有多少个数字等于x。

```C++
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
const int INF=0x3f3f3f3f;
const int maxn=1e6+10;
const int depth=31;

struct Trie{
    #define ls T[x].ch[0]
    #define rs T[x].ch[1]
    int tot;
    struct Node{
        int siz,ch[2],tag;
    } T[maxn<<5];
    void Init(){tot=0;}
    int NewNode(){memset(&T[++tot],0,sizeof(T[0]));return tot;}

    void pushdown(int x)
    {
        int lz=T[x].tag;
        if(lz&1){swap(ls,rs);T[ls].tag++;}
        T[ls].tag+=lz/2; T[rs].tag+=lz/2;
        T[x].tag=0;
    }

    void Insert(int &rt,int x)
    {
        int o=rt?rt:rt=NewNode(),c;
        for(int i=0;i<depth;++i)
        {
            c=x&1; x>>=1; T[o].siz++;
            if(T[o].tag) pushdown(o);
            if(!T[o].ch[c]) T[o].ch[c]=NewNode();
            o=T[o].ch[c];
        }
    }

    int query(int rt,int x,int y)
    {
        int o=rt;
        for(int k=0;k<y;++k)
        {
            if(T[o].tag) pushdown(o);
            o=T[o].ch[x&1];x>>=1;if(!o) break;
        }
        return T[o].siz;
    }

    void Merge(int x,int y)
    {
        T[x].siz+=T[y].siz;
        if(T[x].tag) pushdown(x);
        if(T[y].tag) pushdown(y);
        for(int i=0;i<2;++i)
        {
            if(T[x].ch[i]&&T[y].ch[i]) Merge(T[x].ch[i],T[y].ch[i]);
            if(!T[x].ch[i]&&T[y].ch[i]) T[x].ch[i]=T[y].ch[i];
        }
    }
} trie;
int n,m,rt[maxn],f[maxn];
int find(int x){return f[x]==x?x:f[x]=find(f[x]);}

int main()
{
    while(~scanf("%d",&n))
    {
        scanf("%d",&m);
        memset(rt,0,sizeof rt);
        trie.Init();
        for(int i=1;i<=n;i++)
        {
            f[i]=i;
            int x;scanf("%d",&x);
            trie.Insert(rt[i],x);
        }
        while(m--)
        {
            int op,x,y,z;
            scanf("%d",&op);
            if(op==1)
            {
                scanf("%d%d",&x,&y);
                x=find(x); y=find(y);
                if(x!=y) trie.Merge(rt[x],rt[y]),f[y]=x;
            }
            if(op==2)
            {
                scanf("%d",&x);
                trie.T[rt[find(x)]].tag++;
            }
            if(op==3)
            {
                scanf("%d%d%d",&x,&y,&z);
                x=find(x);
                printf("%d\n",trie.query(rt[x],z,y));
            }
        }
    }
    return 0;
}
```



### **可持久化Trie树+DFS序**

HDU6191给你一棵有根树，树上每个节点有一个值，每次询问以u为根节点的子树异或上x的最大值。 

```c++
const int maxn = 100000 + 10;  
int n, q, cnt, tot;  
int Rank[maxn];//节点在数组上的位置  
int L[maxn];//dfs进入时间戳  
int R[maxn];//dfs出时间戳  
int v[maxn];//顶点权值  
int T[maxn];//n棵字典树的根节点  
struct node{      
    int sum,l, r;      
    node(){ sum = 0; l = 0; r = 0; } 
}Tree[maxn<<6]; 
int head[maxn],ee = 0;  
struct edge{      
    int v , last; 
}Edge[maxn];  
void add(int u, int v)  
{      
    Edge[ee].v = v;      
    Edge[ee].last = head[u];      
    head[u] = ee++;  
} 
void dfs(int u)  
{      
    L[u] = ++cnt;      
    for(int i = head[u]; i != -1 ; i = Edge[i].last)      
    {          
        int v = Edge[i].v;          
        dfs(v);      
    }     
    R[u] = cnt;  
} 
void Insert(int root, int x, int value)  
{      
    int p = root;      
    for(int i = 30; i >= 0; i--)      
    {          
        int num = (x>>i)&1;          
        if(num == 0)          
        {              
            if(Tree[p].l == 0)              
            {                  
                Tree[p].l = tot++;                  
                p = Tree[p].l;                  
                Tree[p].sum += value;              
            }              
            else  p = Tree[p].l,  Tree[p].sum += value;                 
        }          
        else       
        {             
            if(Tree[p].r == 0)             
            {                  
                Tree[p].r = tot++;                  
                p = Tree[p].r;                  
                Tree[p].sum += value;              
            }              
            else              
            {                  
                p = Tree[p].r;                  
                Tree[p].sum += value;              
            }          
        }      
    }  
}  
int update(int pa, int x)  
{      
    int now = tot++;      
    int p = now;     
    for(int i = 30; i >= 0; i--)     
    {          
        int num = (x>>i)&1;          
        if(num == 0)          
        {              
            Tree[p].r = Tree[pa].r;              
            Tree[p].l = tot++;              
            p = Tree[p].l; pa = Tree[pa].l;              
            Tree[p].sum = Tree[pa].sum + 1;          
        }          
        else         
        {              
            Tree[p].l = Tree[pa].l;              
            Tree[p].r = tot++;              
            p = Tree[p].r;  pa = Tree[pa].r;              
            Tree[p].sum = Tree[pa].sum + 1;          
            }      
        }      
        return now;    
}  
int query(int p1, int p2, int x)  
{      
    int d;      
    int ans = 0;      
    for(int i = 30; i >= 0; i--)      
    {          
        int num = (x>>i)&1;          
        if(num == 0)          
        {              
            if(Tree[p2].r == 0)              
            {                  
                p2 = Tree[p2].l;                  
                p1 = Tree[p1].l;                  
                continue;              
            }              
            d = Tree[Tree[p2].r].sum - Tree[Tree[p1].r].sum;              
            if(d > 0)              
            {                  
                ans += (1<<i);                  
                p1 = Tree[p1].r;                  
                p2 = Tree[p2].r;              
            }              
            else  p1 = Tree[p1].l,p2 = Tree[p2].l;                       
        }          
        else          
        {              
            if(Tree[p2].l == 0)              
            {                  
                p2 = Tree[p2].r; p1 = Tree[p1].r;                  
                continue;              
            }              
            d = Tree[Tree[p2].l].sum - Tree[Tree[p1].l].sum;              
            if(d > 0) { ans += (1<<i); p1 = Tree[p1].l; p2 = Tree[p2].l; }              
            else  p1 = Tree[p1].r,p2 = Tree[p2].r;                       
        }      
    }      
    return ans;  
}  
void init()  
{      
    cnt = 0; tot = 0; dfs(1);      
    for(int i = 1; i <= n; i++)  Rank[L[i]] = v[i];       
}  
int main() 
{               
    while(~scanf("%d%d", &n, &q))      
    {          
        for(int i=0; i <(maxn<<6)-1; i++)Tree[i].l=Tree[i].r=Tree[i].sum=0;                  
        for(int i = 1; i <= n; i++)  scanf("%d", &v[i]);                   
        int f;  ee = 0;          
        memset(head, -1, sizeof(head));          
        for(int i = 1; i < n; i++)  scanf("%d", &f),add(f, i + 1);                 
        init();  T[0] = ++tot;          
        for(int i = 1; i <= n; i++)  Insert(T[0], Rank[i], 1);                 
        for(int i = 1; i <= n; i++)  Insert(T[0], Rank[i], -1);                
        for(int i = 1; i <= n; i++)  T[i] = update(T[i - 1], Rank[i]);           
        int uu, xx;          
        for(int i = 1; i <= q; i++)          
        {              
            scanf("%d%d", &uu, &xx);              
            printf("%d\n", query(T[L[uu] - 1], T[R[uu]], xx));            
        }      
    }      
	return 0;  
}  
```



### **01Trie**

删除01Trie中一个数，求和一个数的xor最大最小值

```C++
const int maxn=1e5+10;
int a[maxn];
struct Trie_01{
	int root,tot,nxt[maxn*32][2],cnt[maxn*32],end[maxn*32];
	int Newnode()
	{
		tot++;
		memset(nxt[tot],0,sizeof(nxt[tot]));
		cnt[tot]=0;
		end[tot]=0;
		return tot;
	}
	void Init()
	{
		tot=0;
		root=Newnode();
	}
	void Insert(int x)
	{
		int p=root;
		cnt[p]++;
		for(int i=31;i>=0;--i)
		{
			int idx=(1&(x>>i));
			if(!nxt[p][idx])
				nxt[p][idx]=Newnode();
			p=nxt[p][idx];
			++cnt[p];//可能会有重复的 
		}
		end[p]=x;
	}
	void Delete(int x)//删除一个数x 
	{
		int p=root;
		--cnt[p];
		for(int i=31;i>=0;i--)
		{
			int idx=(1&(x>>i));
			p=nxt[p][idx];
			--cnt[p];
		}
	}
	int QueryMax(int x)//求xor最大值 
	{
		int p=root;
		for(int i=31;i>=0;--i)
		{
			int idx=(1&(x>>i));
			if(idx==0)
			{
				if(nxt[p][1]&&cnt[nxt[p][1]]) p=nxt[p][1]; 
				else p=nxt[p][0];
			}
			else
			{
				if(nxt[p][0]&&cnt[nxt[p][0]]) p=nxt[p][0];
				else p=nxt[p][1];
			}
		}
		return (x^end[p]);
	}
	int QueryMin(int x)//求xor最小值
	{
		int p=root;	
		for(int i=31;i;--i)
		{
			int idx=(1&(x>>i));
			if(idx==1)
			{
				if(nxt[p][1]&&cnt[nxt[p][1]]) p=nxt[p][1]; 
				else p=nxt[p][0];
			}
			else
			{
				if(nxt[p][0]&&cnt[nxt[p][0]]) p=nxt[p][0];
				else p=nxt[p][1];
			}
		}
		return (x^end[p]);
	}
}trie;
```

输入两个数组，可以将这两个数组顺序任意变化，问相同位置异或能生成的字典序最小序列

```C++
const int maxn=1e5+10;
int T,n,a;
vector<pii> ans;
struct Trie_01{
	int root,tot,nxt[maxn*32][2],cnt[maxn*32],end[maxn*32];
	int Newnode()
	{
		++tot;
		memset(nxt[tot],0,sizeof(nxt[tot]));
		cnt[tot]=0;
		end[tot]=0;
		return tot;
	}
	void Init()
	{
		tot=0;
		root=Newnode();
	}
	void Insert(int x)
	{
		int p=root;
		cnt[p]++;
		for(int i=31;i>=0;--i)
		{
			int idx=(1&(x>>i));
			if(!nxt[p][idx])
				nxt[p][idx]=Newnode();
			p=nxt[p][idx];
			++cnt[p];
		}
		end[p]=1;
	}
}trie1,trie2;
void dfs(int x,int y,int num)
{
	int dnum=min(trie1.cnt[x],trie2.cnt[y]);
	trie1.cnt[x]-=dnum;trie2.cnt[y]-=dnum;
	if(trie1.end[x]) {ans.push_back(mkp(num,dnum)); return ;}
	if((trie1.cnt[trie1.nxt[x][0]]) && (trie2.cnt[trie2.nxt[y][0]]))
        dfs(trie1.nxt[x][0], trie2.nxt[y][0], num<<1);
    if((trie1.cnt[trie1.nxt[x][1]]) && (trie2.cnt[trie2.nxt[y][1]]))
        dfs(trie1.nxt[x][1], trie2.nxt[y][1], num<<1);
    if((trie1.cnt[trie1.nxt[x][1]]) && (trie2.cnt[trie2.nxt[y][0]]))
        dfs(trie1.nxt[x][1], trie2.nxt[y][0], (num<<1)+1);
    if((trie1.cnt[trie1.nxt[x][0]]) && (trie2.cnt[trie2.nxt[y][1]]))
        dfs(trie1.nxt[x][0], trie2.nxt[y][1], (num<<1)+1);	
}
int main()
{
	T=read();
	while(T--)
	{
 		trie1.Init();trie2.Init();
 		ans.clear();
		n=read();
		for(int i=1;i<=n;++i) a=read(),trie1.Insert(a);
		for(int i=1;i<=n;++i) a=read(),trie2.Insert(a);
		dfs(trie1.root,trie2.root,0);
		sort(ans.begin(), ans.end());
        for(int i=0,len=ans.size();i<len;++i) 
		{
            for(int j=0,len2=ans[i].se;j<len2;++j) 
			{
                printf("%d",ans[i].fi);
                if(j<ans[i].se-1) printf(" ");
            }
            i<len-1?printf(" "):puts("");
        }
	}
	
	return 0;	
}
```



### **Manacher**

```c++
Len[i]表示以字符str[i](加入符号后的)为中心的最长回文字串的最右字符到T[i]的长度
char str[6000010],s[3000005];  
int n,mx,id,len,Len[3000005];  
void init()
{  
    int k=0; str[k++]='$';  
    for(int i=0;i<len;i++)  str[k++]='#',str[k++]=s[i];  
    str[k++]='#'; len=k;  
}  
int Manacher()
{  
  	Len[0]=0;  
  	int sum=0;mx=0;  
  	for(int i=1;i<len;i++)
	{  
    	if(i<mx) Len[i]=min(mx-i,Len[2*id-i]);  
    	else Len[i]=1;  
   		while(str[i-Len[i]]==str[i+Len[i]]) Len[i]++;  
    	if(Len[i]+i>mx)
		{  
      		mx=Len[i]+i;  id=i;  
      		sum=Max(sum,Len[i]);  
    	}  
  	}  
  	return (sum-1);  
}  
```



### **字符串哈希**

```C++
const int maxn=1010;
const long long mod=1e9+7;
typedef long long ll;
ll p[maxn],h[maxn+1];
set<ll> s; 
int N,base=31;//base如果全部字符base=255,如果全是小数base=131 
char str[maxn];
void pre_work()
{
	p[0]=1;
    for(int i=1;i<maxn;i++) p[i]=p[i-1]*base%mod;
}
void Hash(int n)
{
    h[0]=0;
    for(int i=1;i<=n;i++) h[i]=(h[i-1]*base%mod+str[i]-'a'+1)%mod; 
}
ll get_hash(int l,int r){
	return (h[r]-h[l-1]*p[r-l+1]%mod)%mod;
}
int main()
{
	pre_work();
	scanf("%d",&N);
	while(N--)
	{
		scanf("%s",str+1);
		Hash(strlen(str+1));
		s.insert(get_hash(1,strlen(str+1)));
	}
	printf("%d\n",s.size());
	return 0;
}

```



#### **2019上海网络赛G题**

题意：给了一个母串S, 每次循环给了一个模板串，问模板串在母 串中“匹配”了多少次？“匹配”的意思就是首字母和尾字母一样， 中间字母顺序可以换。

```C++
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef unsigned long long ull;
#define pb push_back
#define base 1000039
#define pii pair<int,int>
#define pil pair<int,ll>
#define mkp make_pair
#define RI register int
const int INF=0x3f3f3f3f;
const int maxn=1e5+10;
const int maxm=2e4+10;
int T,Q,ans[maxm];
char s1[maxn<<1],s2[maxn];
ull p[30],x[maxm],y[maxn];
vector<int> vec[maxn];
inline int idx(char ch){return ch-'a'+1;}
inline void preHash()
{
	p[0]=1;
    for(RI i=1;i<30;++i) 
		p[i]=p[i-1]*base;
}
inline void work()
{
    for(RI i=1;i<=Q;++i)
	{
        scanf("%s",s2+1);
        int len=strlen(s2+1);
        x[i]=idx(s2[1])*p[28]+idx(s2[len])*p[29];
        for(RI j=2;j<=len-1;++j) x[i]+=p[idx(s2[j])];
        vec[len].pb(i);
    }	
}
inline void solve()
{
	int len=strlen(s1+1);
	for(RI i=1;i<=len;++i)
	{
		int siz=vec[i].size();
        if(!siz) continue;
		ull pre=0; int cnt=0;
		for(RI j=1;j<i;++j) pre=pre+p[idx(s1[j])];
        for(RI j=i;j<=len;++j)
		{
            pre=pre+p[idx(s1[j])];
            if(j!=i) pre-=p[idx(s1[j-i])];
            ull res=pre-p[idx(s1[j])]-p[idx(s1[j-i+1])]+(idx(s1[j-i+1]))*p[28]+(idx(s1[j]))*p[29];
            y[++cnt]=res;
        }
        sort(y+1,y+1+cnt);
        for(RI j=0,siz=vec[i].size();j<siz;++j)
		{
            int dn=lower_bound(y+1,y+1+cnt,x[vec[i][j]])-y;
            int up=upper_bound(y+1,y+1+cnt,x[vec[i][j]])-y;
            ans[vec[i][j]]=up-dn;
        }
    }	
}
int main()
{
	preHash(); 
	scanf("%d",&T); 
    while(T--)
	{
		for(RI i=0;i<maxn;++i) vec[i].clear();
        scanf("%s%d",s1+1,&Q);
        work(); solve();
        for(RI i=1;i<=Q;++i) printf("%d\n",ans[i]);
    }
    return 0;
}
```



### **SA(后缀数组)**

#### **最大不重叠相似子串**

```C++
const int INF=0x3f3f3f3f;
const int maxn=20020;
int n,s[maxn];
int rk[maxn],sa[maxn],height[maxn];
int x[maxn<<1],y[maxn<<1],c[maxn];

inline void get_SA(int m) 
{
	for(int i=1;i<=m;++i) c[i]=0;
    for(int i=1;i<=n;++i) ++c[x[i]=s[i]];
    for(int i=2;i<=m;++i) c[i]+=c[i-1];
    for(int i=n;i>=1;--i) sa[c[x[i]]--]=i;
    for(int k=1;k<=n;k<<=1) 
	{
        int num=0;
        for(int i=n-k+1;i<=n;++i) y[++num]=i;
        for(int i=1;i<=n;++i) if(sa[i]>k) y[++num]=sa[i]-k;
	    for(int i=1;i<=m;++i) c[i]=0;
        for(int i=1;i<=n;++i) ++c[x[i]];
        for(int i=2;i<=m;++i) c[i]+=c[i-1]; 
        for(int i=n;i>=1;--i) sa[c[x[y[i]]]--]=y[i],y[i]=0;
        swap(x,y);
        x[sa[1]]=1;
        num=1;
        for(int i=2;i<=n;++i)
            x[sa[i]]=(y[sa[i]]==y[sa[i-1]]&&y[sa[i]+k]==y[sa[i-1]+k])?num:++num;
        if(num==n) break;
        m=num;
    }
}
inline void get_height() 
{
    int k=0;
    for(int i=1;i<=n;++i) rk[sa[i]]=i;
    for(int i=1;i<=n;++i) 
	{
        if(rk[i]==1) continue;
        if(k) --k;
        int j=sa[rk[i]-1];
        while(j+k<=n&&i+k<=n&&s[i+k]==s[j+k]) ++k;
        height[rk[i]]=k;
    }
}
bool check(int k)
{
	int mx=-INF,mi=INF;
	for(int i=1;i<=n;++i)
	{
		if(height[i]>=k)
		{
			mx=max(mx,max(sa[i],sa[i-1]));
			mi=min(mi,min(sa[i],sa[i-1]));
			if(mx-mi>k) return true;
		}	
		else mx=-INF,mi=INF;
	}
	return false;
}

int main()
{
	while(scanf("%d",&n) && n)
	{
		int pre,now; n--;
		scanf("%d",&pre); 
		for(int i=1;i<=n;++i) scanf("%d",&now),s[i]=now-pre+88,pre=now;
		get_SA(176);
		get_height();
		int l=1,r=n>>1;
		while(l+1<r)
		{
			int mid=l+r>>1;
			if(check(mid)) l=mid;
			else r=mid;
		}
		int ans;
		if(check(r)) ans=r;
		else ans=l;
		printf("%d\n",ans>=4? ans+1:0);
	}

	return 0;	
}
```



#### **求两个字符串长度不小于 k 的公共子串的个数**

//论文题,给定两个字符串A和B,求长度不小于 k 的公共子串的个数（可以相同）

```C++
typedef long long ll;
const int MAXN=300000;
int sa[MAXN+9],he[MAXN+9],ra[MAXN+9],xx[MAXN+9],yy[MAXN+9],buc[MAXN+9],q[MAXN+9][2];
char s[MAXN+9];
int len,m;
void get_suf()
{
    int *x=xx,*y=yy;
    for(int i=0;i<m;i++) buc[i]=0;
    for(int i=0;i<len;i++) buc[x[i]=s[i]]++;
    for(int i=1;i<m;i++) buc[i]+=buc[i-1];
    for(int i=len-1;i>=0;i--) sa[--buc[x[i]]]=i;
    for(int k=1;k<=len;k<<=1)
	{
        int p=0;
        for(int i=len-1;i>=len-k;i--) y[p++]=i;
        for(int i=0;i<len;i++) if(sa[i]>=k) y[p++]=sa[i]-k;
        for(int i=0;i<m;i++) buc[i]=0;
        for(int i=0;i<len;i++) buc[x[y[i]]]++;
        for(int i=1;i<m;i++) buc[i]+=buc[i-1];
        for(int i=len-1;i>=0;i--) sa[--buc[x[y[i]]]]=y[i];
        swap(x,y);
        p=1;x[sa[0]]=0;
        for(int i=1;i<len;i++)
		{
            if(y[sa[i-1]]==y[sa[i]]&&y[sa[i-1]+k]==y[sa[i]+k])
                x[sa[i]]=p-1;
            else x[sa[i]]=p++;
        }
        if(p>=len) break;
        m=p;
    }
    for(int i=0;i<len;i++) ra[sa[i]]=i;
    int k=0;
    for(int i=0;i<len;i++)
	{
        if(ra[i]==0) { he[0]=0; continue; }
        if(k) k--;
        int j=sa[ra[i]-1];
        while(s[i+k]==s[j+k]&&i+k<len&&j+k<len) k++;
        he[ra[i]]=k;
    }
}
ll solve(int len1,int k)
{
    ll ans=0,cnt=0,sum=0,top=0;
    for(int i=1;i<len;i++)
	{
        if(he[i]<k) { top=sum=0;continue; }
        cnt=0;
        if(sa[i-1]<len1) { cnt++;sum+=he[i]-k+1; }
        while(top&&he[i]<=q[top][1])
		{
            sum-=q[top][0]*(q[top][1]-he[i]);
            cnt+=q[top--][0];
        }
        q[++top][0]=cnt;
        q[top][1]=he[i];
        if(sa[i]>len1) ans+=sum;
    }
    sum=0;top=0;
    for(int i=1;i<len;i++)
	{
        if(he[i]<k) { top=sum=0;continue; }
        cnt=0;
        if(sa[i-1]>len1) { cnt++;sum+=he[i]-k+1; }
        while(top&&he[i]<=q[top][1])
		{
            sum-=q[top][0]*(q[top][1]-he[i]);
            cnt+=q[top--][0];
        }
        q[++top][0]=cnt;
        q[top][1]=he[i];
        if(sa[i]<len1) ans+=sum;
    }
    return ans;
}
int main()
{
    int k;
    while(scanf("%d",&k)&&k)
	{
        scanf("%s",s);
        int len1=strlen(s);
        s[len1]='#';
        scanf("%s",s+len1+1);
        len=strlen(s);
        m=200;
        get_suf();
        printf("%lld\n",solve(len1,k));
    }
    return 0;
}
```



### **SAM (后缀自动机)**

后缀自动机的性质

1.有一个源点，边代表在当前字符串后增加一个字符。

​	2.每个点代表一个 endpos 等价类，到达一个点的路径形成的子串必须属于此点的类。

​	3.点之间有父子关系，到达点 i 的所有字符串的长度都必然大于到达 fa(i) 的所有字符串的长度，且到达 fa(i) 的任意一字符串必为到达 i 的任意一字符串的后缀。

Sam 讲解：<https://www.luogu.org/blog/Kesdiael3/hou-zhui-zi-dong-ji-yang-xie>

parent tree上DP得到每一个类中的子串**和在原串中出现的位置相关**的一些信息。例如，出现的次数，第一次出现的位置，出现在某个位置之后的次数等。

#### **洛谷p3975 求字典序第K小串** 

```c++
struct SAM{//求字典序第K小串 
	int l[maxn<<1],fa[maxn<<1],nxt[maxn<<1][26];
	int last,cnt,c[maxn<<1],siz[maxn<<1],sum[maxn<<1],a[maxn<<1];
	void Init()
	{
		memset(siz,0,sizeof(siz));
		memset(c,0,sizeof(c));
		memset(sum,0,sizeof(sum));
		memset(a,0,sizeof(a));
		last=cnt=1;
		memset(nxt[1],0,sizeof(nxt[1]));
		fa[1]=l[1]=0;
	}
	int NewNode()
	{
		cnt++;
		memset(nxt[cnt],0,sizeof(nxt[cnt]));
		fa[cnt]=l[cnt]=0;
		return cnt;
	}
	void Add(int ch)
	{
		int p=last,np=NewNode();
		last=np; l[np]=l[p]+1;
		siz[np]=1;
		while(p&&!nxt[p][ch]) nxt[p][ch]=np,p=fa[p];
		if(!p) fa[np]=1;
		else
		{
			int q=nxt[p][ch];
			if(l[q]==l[p]+1) fa[np]=q;
			else
			{
				int nq=NewNode();
				memcpy(nxt[nq],nxt[q],sizeof(nxt[q]));
				fa[nq]=fa[q];
				l[nq]=l[p]+1;
				fa[np]=fa[q]=nq;
				while(nxt[p][ch]==q) nxt[p][ch]=nq,p=fa[p];
			}
		}
	}
	void Build()
	{
        int len=strlen(s+1);
        for(int i=1;i<=len;i++) Add(s[i]-'a');
	}
//拓扑排序实现 长度(longest)对应的状态i(endpos)
//a[x]:长度第x小的子串对应的状态
	void topusort()
	{
		for(int i=1;i<=cnt;++i) c[l[i]]++;
		for(int i=1;i<=cnt;++i) c[i]+=c[i-1];
		for(int i=1;i<=cnt;++i) a[c[l[i]]--]=i;
		for(int i=cnt;i;--i)
		{//t==1:不同位置的相同子串视为不同 
			if(t) siz[fa[a[i]]]+=siz[a[i]];
			else siz[a[i]]=1;//t==0:视为相同 
		}
		siz[1]=0;
		for(int i=cnt;i;--i)
		{
			sum[a[i]]=siz[a[i]];
			for(int j=0;j<26;++j)
				if(nxt[a[i]][j]) sum[a[i]]+=sum[nxt[a[i]][j]];
		}
	}	
	void dfs()
	{
		if(k>sum[1]){puts("-1");return ;}
		int now=1;
		while(k>0)
		{
			int p=0;
			while(k>sum[nxt[now][p]])
			{
				k-=sum[nxt[now][p]];
				p++;
			}
			now=nxt[now][p];
			putchar('a'+p);
			k-=siz[now];
		}
		return ;
	}
} sam;
int main()
{
	scanf("%s%d%d",s+1,&t,&k);
	sam.Init();   sam.Build(); 
	sam.topusort();  sam.dfs();
	return 0;
}
```



#### **动态求出现至少k次本质不同子串个数**

```c++
#include<bits/stdc++.h>  
using namespace std;  
const int maxn=2e5+1000;  
char s[maxn];  
int len,k,n,m;  
char temp[5];  
struct SAM{  
    int last,cnt,nxt[maxn<<1][26],fa[maxn<<1],l[maxn<<1],num[maxn<<1];  
    int ans;  
    void init()
	{  
        last=cnt=1;  
        memset(nxt[1],0,sizeof nxt[1]);  
        fa[1]=l[1]=num[1]=0;  
        ans=0;  
    }  
    int newnode()
	{  
        cnt++;  
        memset(nxt[cnt],0,sizeof nxt[cnt]);  
        fa[cnt]=l[cnt]=num[cnt]=0;  
        return cnt;  
    }  
    void add(int c)
	{  
        int p=last,np=newnode();    
        last=np;l[np] =l[p]+1;    
        while(p&&!nxt[p][c]) nxt[p][c]=np,p=fa[p];       
        if(!p) fa[np]=1;
		else
		{    
            int q=nxt[p][c];    
            if(l[q]==l[p]+1) fa[np]=q;
			else
			{    
                int nq=newnode();    
                memcpy(nxt[nq],nxt[q],sizeof nxt[q]);    
                fa[nq]=fa[q]; 
				num[nq]=num[q];////  
                l[nq]=l[p]+1; 
				fa[np]=fa[q]=nq;    
                while(nxt[p][c]==q) nxt[p][c]=nq,p=fa[p];       
            }    
        }  
        int temp=last;  
        while(temp)
		{  
            if(num[temp]>=k) break; 
            num[temp]++;  
            if(num[temp]==k) ans+=l[temp]-l[fa[temp]];  
            temp=fa[temp];  
        }    
    }  
} sam;  
int main()
{  
    while(scanf("%d%d%d",&n,&m,&k)!=EOF) 
	{  
        scanf("%s",s);  
        len=strlen(s);  
		sam.init();  
        for(int i=0;i<len;i++) sam.add(s[i]-'a');  
        while(m--)
		{  
            int flag;  
            scanf("%d",&flag);  
            if(flag==1)
			{  
                scanf("%s",temp);  
                sam.add(temp[0]-'a');  
            }
			else printf("%d\n",sam.ans);   
        }  
    }  
    return 0;  
}  
```



#### **线段树合并:求在串s的l,r区间的子串第k个出现位置** 

```c++
const int maxn=1e5+10;
char s[maxn];
int T[maxn<<1][30],fa[maxn<<1],len[maxn<<1],cnt,last,n;
void Insert(int v)
{
    int i,p=last,np,q,nq;
    last=np=++cnt;
    len[np]=len[p]+1;
    for(;p&&T[p][v]==0;p=fa[p]) T[p][v]=np;
    if(p==0) fa[np]=1;
    else
	{
        q=T[p][v];
        if(len[q]==len[p]+1) fa[np]=q;
        else
		{
            nq=++cnt;
            len[nq]=len[p]+1;
            fa[nq]=fa[q];
            fa[q]=fa[np]=nq;
            for(i=0;i<27;i++) T[nq][i]=T[q][i];
            for(;T[p][v]==q;p=fa[p]) T[p][v]=nq;
        }
    }
}
int root[maxn<<1],num;
struct Q{
    int L,R,sum;
}A[maxn*40];
void update(int &x,int l,int r,int k,int v)
{
    A[++num]=A[x];
    x=num;
    A[x].sum++;
    if(l==r) return;
    int mid=(l+r)/2;
    if(k<=mid) update(A[x].L,l,mid,k,v);
    else update(A[x].R,mid+1,r,k,v);
}
int mer(int a,int b,int l,int r)
{
    if(a==0||b==0) return a+b;
    int z=++num,mid=(l+r)/2;
    if(l==r)
	{
         A[z].sum=A[a].sum|A[b].sum;
         return z;
    }
    A[z].L=mer(A[a].L,A[b].L,l,mid);
    A[z].R=mer(A[a].R,A[b].R,mid+1,r);
    A[z].sum=A[A[z].L].sum+A[A[z].R].sum;
    return z;
}
int qkth(int x,int l,int r,int k)
{
    if(l==r) return l;
    int mid=(l+r)/2;
    if(k<=A[A[x].L].sum) return qkth(A[x].L,l,mid,k);
    else return qkth(A[x].R,mid+1,r,k-A[A[x].L].sum);
}
vector<int>g[maxn<<1];
int Fa[maxn<<1][25],pos[maxn<<1];
void dfs(int u)
{
    int i,v;
    for(i=1;i<=19;i++) Fa[u][i]=Fa[Fa[u][i-1]][i-1];
    for(i=0;i<g[u].size();i++)
	{
        v=g[u][i];
        Fa[v][0]=u;
        dfs(v);
        root[u]=mer(root[u],root[v],1,n);
    }
}
int main()
{
    int i,m,t;
    scanf("%d",&t);
    while(t--)
	{
        scanf("%d%d",&n,&m);
        scanf("%s",s+1);
        cnt=last=1;
        memset(T,0,sizeof(T));
        for(i=0;i<=200000;i++) g[i].clear();

        for(i=1;s[i];i++) Insert(s[i]-'a');
        for(i=1;i<=cnt;i++) g[fa[i]].push_back(i);
        int p=1,v;

        memset(root,0,sizeof(root));
        A[0].L=A[0].R=A[0].sum=num=0;
        for(i=1;s[i];i++)
		{
            v=s[i]-'a';
            p=T[p][v];
            pos[i]=p;
            update(root[p],1,n,i,1);
        }
        dfs(1);//合并子树 
        int l,r,k,u,a;
        while(m--)
		{
            scanf("%d%d%d",&l,&r,&a);
            p=pos[r];
            k=r-l+1;
            for(i=19;i>=0;i--) if(len[Fa[p][i]]>=k) p=Fa[p][i];
            if(a>A[root[p]].sum) printf("-1\n");
            else
			{
                u=qkth(root[p],1,n,a);
                printf("%d\n",u-k+1);
            }
        }
    }
    return 0;
}
```



#### **两个字符串的子串拼接成的不同字符串数量**

给你两个字符串，然后从第一个字符串里面取出一个子串X，从第二个字符串里面取出一个子串Y，两个拼接在一起组成新的字符串，其中X、Y都可以是空串，问有多少个这样不同的串。

```c++
typedef unsigned long long ull;
#define RI register int
const int maxn=1e5+10;
char s1[maxn],s2[maxn];
struct SAM{
	int last,tot,nxt[maxn<<1][27],fa[maxn<<1],l[maxn<<1];
	inline void Init()
	{
		last=tot=1;
		memset(nxt[tot],0,sizeof(nxt[tot]));
		l[tot]=fa[tot]=0;
	}
	inline int NewNode()
	{
		++tot;
		memset(nxt[tot],0,sizeof(nxt[tot]));
		l[tot]=fa[tot]=0;
		return tot;
	}
	inline void Add(int c)
	{
		int np=NewNode(),p=last;
		last=np;l[np]=l[p]+1;
		while(p&&!nxt[p][c]) nxt[p][c]=np,p=fa[p];
		if(!p) fa[np]=1;
		else
		{
			int q=nxt[p][c];
			if(l[q]==l[p]+1) fa[np]=q;
			else
			{
				int nq=NewNode();
				memcpy(nxt[nq],nxt[q],sizeof(nxt[q]));
				fa[nq]=fa[q];
				l[nq]=l[p]+1;
				fa[q]=fa[np]=nq;
				while(p&&nxt[p][c]==q) nxt[p][c]=nq,p=fa[p];
			}
		}
	}
} sam1,sam2;

int T;
ull dp1[maxn<<1],dp2[maxn<<1];
inline ull dfs2(int u)
{
	if(!u) return 0;
	if(dp2[u]) return dp2[u];
	ull res=1;
	for(int i=0;i<26;++i)
	{
		int nt=sam2.nxt[u][i];
		if(nt) res+=dfs2(nt);	
	}
	return dp2[u]=res;	
}
inline ull dfs(int u)
{
	if(dp1[u]) return dp1[u];
	ull res=1;
	for(int i=0;i<26;++i)
	{
		int nt=sam1.nxt[u][i];
		if(nt) res+=dfs(nt);
		else res+=dfs2(sam2.nxt[1][i]);
	}
	return dp1[u]=res;
}

int main()
{
	scanf("%d",&T);
	while(T--)
	{
		sam1.Init();sam2.Init();
		memset(dp1,0,sizeof dp1);
		memset(dp2,0,sizeof dp2);
		
		scanf("%s%s",s1,s2);
		for(int i=0;s1[i];++i) sam1.Add(s1[i]-'a');
		for(int i=0;s2[i];++i) sam2.Add(s2[i]-'a');
		
		printf("%I64u\n",dfs(1));
	}

	return 0;	
}
```



#### **SAM+线性基**

HDU6694给你一个字符串S，然后q个询问，每次给出S的一个子串T。对于每个询问的子串T，Calabash可以在Ｓ中选择任意个以T作为后缀的子串，然后生成子串对应数目个石子堆，每堆的石子数量等于w[对应子串在S中出现的次数]。然后Rounddog可以从这么多堆石子中选择任意堆的石子（至少选一堆），两人开始玩Nim游戏，Calabash先手。现在问Calabash是否存在必胜策略，如果有输出Calabash在必胜策略下，能够选出来的最多的石子数目。

 

```c++
#include<bits/stdc++.h>
#define INF 0x3f3f3f3f
#define eps 1e-4
#define pi 3.141592653589793
#define P 1000000007
#define LL long long
#define ULL unsigned LL
#define pb push_back
using namespace std;
const int N = 200010;
int n,f[N],num[N],dep[N],dp[19][N],Right[N],tot;
vector<int>g[N]; LL w[N],p[N];
char s[N];

struct Linear_Basis{
    LL b[58]; ULL ans;
    inline void init(){ans=0;memset(b,0,sizeof(b));}
    inline bool ins(LL x)
    {
        LL tmp=x;
        for(int i=57;i>=0;i--)
            if (x&(1LL<<i))
            {
                if (!b[i]) {b[i]=x;ans+=tmp;break;}
                x^=b[i];
            }
        return x>0;
    }
} LB[N];
 
struct Suffix_Automation{
    int tot,cur;
    struct node{int ch[26],len,fa;} T[N];
    inline void init()
    {
        cur=tot=1;memset(T,0,sizeof(T));
        memset(Right,0,sizeof(Right));
    }
    void ins(int x,int id)
    {
        int p=cur;cur=++tot;
		T[cur].len=id;Right[cur]++;
        for(;p&&!T[p].ch[x];p=T[p].fa) T[p].ch[x]=cur;
        if (!p) {T[cur].fa=1;return;}
		int q=T[p].ch[x];
        if (T[p].len+1==T[q].len) {T[cur].fa=q;return;}
        int np=++tot; memcpy(T[np].ch,T[q].ch,sizeof(T[q].ch));
        T[np].fa=T[q].fa; T[q].fa=T[cur].fa=np; T[np].len=T[p].len+1;
        for(;p&&T[p].ch[x]==q;p=T[p].fa) T[p].ch[x]=np;
    }
    inline void build()
    {
        for(int i=2;i<=tot;i++) g[T[i].fa].pb(i);
    }
 
} SAM;
 
void ST()
{
    for(int j=1;j<19;j++)
        for(int i=1;i<=tot;i++)
            dp[j][i]=dp[j-1][dp[j-1][i]];
}
 
void dfs(int x,int fa,int depth)
{
    dp[0][x]=fa;
    dep[x]=depth;
    for (auto i:g[x])
    {
        if (i==fa) continue;
        dfs(i,x,depth+1);
        Right[x]+=Right[i];
    }
}
 
namespace IO{
    #define BUF_SIZE 100000
    #define OUT_SIZE 100000
    #define ll long long
    //fread->read
    bool IOerror=0;
    inline char nc(){
        static char buf[BUF_SIZE],*p1=buf+BUF_SIZE,*pend=buf+BUF_SIZE;
        if (p1==pend){
            p1=buf; pend=buf+fread(buf,1,BUF_SIZE,stdin);
            if (pend==p1){IOerror=1;return -1;}
        }
        return *p1++;
    }
    inline bool blank(char ch){return ch==' '||ch=='\n'||ch=='\r'||ch=='\t';}
    inline void read(int &x){
        bool sign=0; char ch=nc(); x=0;
        for (;blank(ch);ch=nc());
        if (IOerror)return;
        if (ch=='-')sign=1,ch=nc();
        for (;ch>='0'&&ch<='9';ch=nc())x=x*10+ch-'0';
        if (sign)x=-x;
    }
    inline void read(ll &x){
        bool sign=0; char ch=nc(); x=0;
        for (;blank(ch);ch=nc());
        if (IOerror)return;
        if (ch=='-')sign=1,ch=nc();
        for (;ch>='0'&&ch<='9';ch=nc())x=x*10+ch-'0';
        if (sign)x=-x;
    }
    inline void read(char *s){
        char ch=nc();
        for (;blank(ch);ch=nc());
        if (IOerror)return;
        for (;!blank(ch)&&!IOerror;ch=nc())*s++=ch;
        *s=0;
    }
}
 
struct Ostream_fwrite{
        char *buf,*p1,*pend;
        Ostream_fwrite(){buf=new char[BUF_SIZE];p1=buf;pend=buf+BUF_SIZE;}
        void out(char ch){
            if (p1==pend){
                fwrite(buf,1,BUF_SIZE,stdout);p1=buf;
            }
            *p1++=ch;
        }
        void print(ULL x){
            static char s[20],*s1;s1=s;
            if (!x)*s1++='0';if (x<0)out('-'),x=-x;
            while(x)*s1++=x%10+'0',x/=10;
            while(s1--!=s)out(*s1);
        }
        void flush(){if (p1!=buf){fwrite(buf,1,p1-buf,stdout);p1=buf;}}
        void print(char *s){while (*s)out(*s++);}
        ~Ostream_fwrite(){flush();}
}Ostream;
inline void print(ULL x){Ostream.print(x);}
inline void print(char *s){Ostream.print(s);}
inline void flush(){Ostream.flush();}
inline bool cmp(int a,int b){return p[a]>p[b];}
using namespace IO;
 
 
int main()
{
    int T; read(T);
    while(T--)
    {
        SAM.init();
        memset(g,0,sizeof(g));
        read(n); read(s);
        for(int i=0;i<n;i++)
        {
            SAM.ins(s[i]-'a',i+1);
            num[i+1]=SAM.cur;
        }
        tot=SAM.tot;
        for(int i=1;i<=n;i++) read(w[i]);
        SAM.build(); dfs(1,0,0); ST();
        for(int i=1;i<=tot;i++)
            LB[i].init(),f[i]=i,p[i]=w[Right[i]];
        sort(f+1,f+1+tot,cmp);
        for(int i=1;i<=tot;i++)
            for(int now=f[i];now;now=dp[0][now])
                if (!LB[now].ins(p[f[i]])) break;
        int q; read(q);
        while(q--)
        {
            int x,y;
            read(x); read(y);
            int len=y-x+1;y=num[y];
            for(int i=18;i>=0;i--)
                if (SAM.T[dp[i][y]].len>=len) y=dp[i][y];
            print(LB[y].ans); print("\n");
        }
        flush();
    }
    return 0;
}
```



### **GSAM(广义后缀自动机)**

#### **一颗字典树每次查询一个给出字符串是字典树上多少串的后缀**

```c++
//build sam using trie
#include<bits/stdc++.h>
using namespace std;
const int maxn = 1e6+100;
typedef long long ll;
struct Suffix_Automaton{
    int nxt[maxn*2][26],fa[maxn*2],l[maxn*2];
    int last,cnt;
    vector<int> E[maxn*2];
    int Num[maxn*2];
    Suffix_Automaton(){ clear(); }
    void clear()
	{
        last =cnt=1;
        fa[1]=l[1]=0;
        memset(nxt[1],0,sizeof nxt[1]);
    }
    int add(int pre,int c,int num)
	{
        last = pre;
        int p = last;
        int np = ++cnt;
        Num[np] = num;
        memset(nxt[cnt],0,sizeof nxt[cnt]);
        l[np] = l[p]+1;last = np;
        while (p&&!nxt[p][c])nxt[p][c] = np,p = fa[p];
        if(!p) fa[np]=1;
        else
		{
            int q = nxt[p][c];
            if (l[q]==l[p]+1)fa[np] =q;
            else
			{
                int nq = ++ cnt;
                l[nq] = l[p]+1;
                memcpy(nxt[nq],nxt[q],sizeof (nxt[q]));
                fa[nq] =fa[q];fa[np] = fa[q] =nq;
                while (nxt[p][c]==q)nxt[p][c] =nq,p = fa[p];
            }
        }
        return np;
    }
    int dfsl[maxn*2],dfsr[maxn*2];
    int dfn = 0;
    ll sum[maxn*2];
    void dfs(int u)
	{
        dfsl[u]=++dfn;
        sum[dfn]=Num[u];
        for(int v : E[u]) dfs(v);
        dfsr[u]=dfn;
    }
    void build()
	{
        for(int i=2;i<=cnt;i++) E[fa[i]].push_back(i);
        dfs(1);
        for(int i=1;i<=cnt;i++) sum[i]+=sum[i-1];
    }
    void query(char * s)
	{
        int temp = 1;
        while(*s)
		{
            int ch=*s-'A';
            if(!nxt[temp][ch])
			{
                printf("0\n");
                return;
            }
            temp=nxt[temp][ch];
            s++;
        }
        ll ans=sum[dfsr[temp]]-sum[dfsl[temp]-1];
        printf("%lld\n",ans);
    }
}sam;
struct Trie{
    int Root=1;
    int cnt=2;
    int nxt[maxn][26];
    int num[maxn];
    int sam_pos[maxn];
    int add(int p,int ch)
	{
        if(!nxt[p][ch]) nxt[p][ch]=cnt++;
        int now=nxt[p][ch];
        num[now]++;
        return now;
    }
    void bfs()
	{
        queue<int> Q;
        Q.push(1);
        sam_pos[1]=1;
        while(!Q.empty())
		{
            int head=Q.front();
            Q.pop();
            for(int i=0;i<26;i++)
			{
                if(!nxt[head][i])continue;
                int now=nxt[head][i];
                sam_pos[now]=sam.add(sam_pos[head],i,num[now]);
                Q.push(now);
            }
        }
    }
}trie;
int trie_pos[maxn];
int main()
{
    int n,k;
    scanf("%d%d",&n,&k);
    trie_pos[0]=1;
    for(int i=1;i<=n;i++)
	{
        static char s[5];
        int p;
        scanf("%s%d",s,&p);
        int ch=s[0]-'A';
        trie_pos[i]=trie.add(trie_pos[p],ch);
    }
    trie.bfs();
    sam.build();
    for(int i=0;i<k;i++)
	{
        static char t[maxn];
        scanf("%s",t);
        int N=strlen(t);
        reverse(t,t+N);
        sam.query(t);
    }
    return 0;
}
```



#### **线段树合并**

给你n个单词，Q个询问: l ,r ,x在第l到r个单词之间，第x个单词出现了多少次。CF547E(CSL)

```c++
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
const int INF=0x3f3f3f3f;
const int maxn=2e5+10;
int N,Q;
char s[maxn];

int nxt[maxn<<1][26],l[maxn<<1],fa[maxn<<1];
int last,tot,cnt[maxn<<1],c[maxn<<1];
int sz,p[maxn],T[maxn<<1];

struct Tr{
    int ls,rs;
    int num;
} tr[maxn*40];

void Update(int &x,int l,int r,int pos)
{
    if(!x) x=++sz;
    tr[x].num++;
    if(l==r) return ;
    int mid=l+r>>1;
    if(pos<=mid) Update(tr[x].ls,l,mid,pos);
    else Update(tr[x].rs,mid+1,r,pos);
}

int Query(int x,int l,int r,int L,int R)
{
    if(!x) return 0;
    if(L<=l&&r<=R) return tr[x].num;
    int mid=l+r>>1,res=0;
    if(L<=mid) res+=Query(tr[x].ls,l,mid,L,R);
    if(R>mid) res+=Query(tr[x].rs,mid+1,r,L,R);
    return res;
}

int Merge(int x,int y)
{
    if(!x||!y) return x+y;
    int z=++sz;
    tr[z].ls=Merge(tr[x].ls,tr[y].ls);
    tr[z].rs=Merge(tr[x].rs,tr[y].rs);
    tr[z].num=tr[x].num+tr[y].num;
    return z;
}

void Mer()
{
    for(int i=1;i<=tot;++i) cnt[l[i]]++;
    for(int i=1;i<=tot;++i) cnt[i]+=cnt[i-1];
    for(int i=1;i<=tot;++i) c[cnt[l[i]]--]=i;
    for(int i=tot,x;i>1;--i) x=c[i],T[fa[x]]=Merge(T[x],T[fa[x]]);
}

void Init()
{
    last=tot=1; sz=0;
    memset(nxt[tot],0,sizeof nxt[tot]);
    l[tot]=fa[tot]=0;
}

int NewNode()
{
    ++tot;
    memset(nxt[tot],0,sizeof nxt[tot]);
    l[tot]=fa[tot]=0;
    return tot;
}

void Insert(int ch,int x)
{
    int p,q,np,nq;
    if(nxt[last][ch])
    {
        p=last;q=nxt[p][ch];
        if(l[q]==l[p]+1) last=q;//////
        else
        {
            nq=NewNode();
            l[nq]=l[p]+1;fa[nq]=fa[q];
            memcpy(nxt[nq],nxt[q],sizeof(nxt[q]));
            fa[q]=nq;
            while(p&&nxt[p][ch]==q) nxt[p][ch]=nq,p=fa[p];
            last=nq;//////
        }
    }
    else
    {
        np=NewNode(),p=last;
        last=np; l[np]=l[p]+1;
        while(p&&!nxt[p][ch]) nxt[p][ch]=np,p=fa[p];
        if(!p) fa[np]=1;
        else
        {
            q=nxt[p][ch];
            if(l[q]==l[p]+1) fa[np]=q;
            else
            {
                nq=NewNode();
                memcpy(nxt[nq],nxt[q],sizeof nxt[q]);
                fa[nq]=fa[q];
                l[nq]=l[p]+1;
                fa[q]=fa[np]=nq;
                while(p&&nxt[p][ch]==q) nxt[p][ch]=nq,p=fa[p];
            }
        }
    }
}


int main()
{
    scanf("%d%d",&N,&Q);
    Init();
    for(int i=1;i<=N;++i)
    {
        scanf("%s",s);
        int len=strlen(s); last=1;
        for(int j=0;j<len;++j) Insert(s[j]-'a',i),Update(T[last],1,N,i);
        p[i]=last;
    }
    Mer();
    while(Q--)
    {
        int l,r,x;
        scanf("%d%d%d",&l,&r,&x);
        printf("%d\n",Query(T[p[x]],1,N,l,r));
    }
    return 0;
}
```



#### **长度<=m的子串的期望**

HDU6405

给你n个字符串，然后每个字符串有一个快乐值。然后给你m个询问，每个询问给你一个长度，让你写出一个不大于这个长度的字串。这个字串的权值定义为，如果这个字符串中出现过第i个给定字符串的子串，那么权值乘以第i个字符串的快乐值，最后答案就是多个快乐值相乘。现在问你给定长度的字符串权值的期望。

```c++
#include<bits/stdc++.h>
using namespace std;
#define mod 1000000007
typedef long long ll;
const int maxn=1e6+10;
int pw[maxn],n,m;
int flag[maxn],h[maxn],sum[maxn],ans[maxn];
string s[maxn];
bool vis[maxn];
int qpow(int a,int b)
{
    int ans=1;
    while(b)
    {
        if(b&1)ans=(ll)ans*a%mod;
        a=(ll)a*a%mod; b>>=1;
    }
    return ans;
}
void Init()
{
    pw[0]=1;
    for(int i=1;i<maxn;i++)
        pw[i]=(ll)pw[i-1]*26ll%mod;
    for(int i=2;i<maxn;i++)
        pw[i]=(pw[i]+pw[i-1])%mod;
}
struct SAM{
	int fa[maxn<<1],l[maxn<<1],nxt[maxn<<1][26];
	int last,tot;
	void Init()
	{
	    last=tot=1; 
	    memset(nxt[tot],0,sizeof nxt[tot]);
	    l[tot]=fa[tot]=0;
	}
	int NewNode()
	{
	    ++tot;
	    memset(nxt[tot],0,sizeof nxt[tot]);
	    l[tot]=fa[tot]=0;
	    return tot;
	}	
	void Insert(int ch)
	{
	    int p,q,np,nq;
	    if(nxt[last][ch])
	    {
	        p=last;q=nxt[p][ch];
	        if(l[q]==l[p]+1) last=q;//////
	        else
	        {
	            nq=NewNode();
	            l[nq]=l[p]+1;fa[nq]=fa[q];
	            memcpy(nxt[nq],nxt[q],sizeof(nxt[q]));
	            fa[q]=nq;
	            while(p&&nxt[p][ch]==q) nxt[p][ch]=nq,p=fa[p];
	            last=nq;//////
	        }
	    }
	    else
	    {
	        np=NewNode(),p=last;
	        last=np; l[np]=l[p]+1;
	        while(p&&!nxt[p][ch]) nxt[p][ch]=np,p=fa[p];
	        if(!p) fa[np]=1;
	        else
	        {
	            q=nxt[p][ch];
	            if(l[q]==l[p]+1) fa[np]=q;
	            else
	            {
	                nq=NewNode();
	                memcpy(nxt[nq],nxt[q],sizeof nxt[q]);
	                fa[nq]=fa[q];
	                l[nq]=l[p]+1;
	                fa[q]=fa[np]=nq;
	                while(p&&nxt[p][ch]==q) nxt[p][ch]=nq,p=fa[p];
	            }
	        }
	    }
	}
	void cal(string s,int val,int tag)
    {
        int cur=1;
        for(int i=0;s[i];i++)
        {
            cur=nxt[cur][s[i]-'a'];
            for(int tmp=cur;tmp&&flag[tmp]!=tag;tmp=fa[tmp])
                ans[tmp]=1LL*ans[tmp]*val%mod,flag[tmp]=tag;
        }
    }
    void build(int cur)
    {
        vis[cur]=1;
        sum[l[fa[cur]]+1]=(sum[l[fa[cur]]+1]+ans[cur])%mod;
        sum[l[cur]+1]=(sum[l[cur]+1]-ans[cur]+mod)%mod;
        for(int i=0;i<26;i++)
        {
            int Nxt=nxt[cur][i];
            if(Nxt&&!vis[Nxt]) build(Nxt);
        }
    }
} sam;
int main()
{
	ios::sync_with_stdio(0);
	cin.tie(0); cout.tie(0);
	Init(); sam.Init();
	cin>>n;
	for(int i=1;i<=n;++i)
	{
		cin>>s[i]; sam.last=1;
		for(int j=0;s[i][j];++j)
			sam.Insert(s[i][j]-'a');
	}
	fill(ans,ans+maxn,1);
	for(int i=1;i<=n;++i) cin>>h[i];
    for(int i=1;i<=n;++i) sam.cal(s[i],h[i],i);
    sam.build(1); sum[0]=0;
	for(int i=1;i<maxn;++i)//sum[i]表示所有长度为i的串的贡献 
		sum[i]=(sum[i]+sum[i-1])%mod;
    for(int i=1;i<maxn;++i)//sum[i]就表示所有长度为1~i的串的贡献
		sum[i]=(sum[i]+sum[i-1])%mod;
    cin>>m;
    while(m--)
    {
        int k; cin>>k;
        cout<<1LL*sum[k]*qpow(pw[k],mod-2)%mod<<endl;
    }
	return 0;	
} 
```

 

### **ACAM (AC自动机)**

#### **HDU2222:查找模式串**

fail含义是指以x为结尾的后缀在其他模式串中所能匹配的最长前缀的长度

​	

```c++
#include<bits/stdc++.h>
	using namespace std;
	#define maxn 500005  
	int n,T;  
	char s[maxn<<1];   								
    int ch[maxn][26],fail[maxn],val[maxn],last[maxn]; 
    int times[10005],cnt[10005]; // 每个单词出现次数   
    int tot,num;                //单词总数以及单词出现了几个  
    void init()  
    {  
        num=0;  tot=1;  
        memset(ch[0],0,sizeof(ch[0]));  
        memset(cnt,0,sizeof(cnt));  
        memset(times,0,sizeof(times));  
    }     
	int idx(char c){ return c - 'a';}  
    void insert(char*s,int v) //插入
    {  
        int u=0,n=strlen(s);  
        for(int i=0;i<n;++i)  
        {  
            int c=idx(s[i]);  
            if(!ch[u][c])  
            {  
                memset(ch[tot],0,sizeof(ch[tot]));  
                val[tot]=0;  
                ch[u][c]=tot++;   
            }  
            u=ch[u][c];  
        }  
        if(val[u]) times[val[u]]++;//重复的模式串  
        else  val[u]=v,times[v]=1;  
    }
	void getFail()  
    {  
        queue<int> q;  
        fail[0]=0;  
        for(int c=0;c<26;++c)  
        {  
            int u=ch[0][c];  
            if(u){ fail[u]=0;q.push(u);last[u]=0; }  
        }  
        while(!q.empty())  
        {  
            int r=q.front(); q.pop();  
            for(int c=0;c<26;++c)  
            {  								
                int u=ch[r][c];  
                if(!u){ch[r][c]=ch[fail[r]][c];continue;}  
                q.push(u);  
                int v=fail[r];  
                fail[u]=ch[v][c];  
                last[u] = val[fail[u]]?fail[u]:last[fail[u]];  
            }  
        }  
    }   
    void print(int i,int j)                   
    {  
        if(j)   
        {  
            if(!cnt[val[j]])  num+=times[val[j]];  
            cnt[val[j]]++; //每个模式串出现的次数 
            print(i,last[j]);  
        }  
    }  
    void find(char *T)  
    {  
        int n=strlen(T);  
        int j=0;  
        for(int i=0;i<n;++i)  
        {  
            int c=idx(T[i]);  
          	j=ch[j][c];  
            if(val[j]) print(i,j);  
            else if(last[j]) print(i,last[j]);  
        }  
    }    
	int main()  
	{  
	    cin>>T;
	    while(T--)  
	    {  
	        init(); cin>>n;  
	        for(int i=1;i<=n;i++)  
	        {  
	            scanf("%s",s);  
	            insert(s,i);  
	        }  
	        getFail();  
	        scanf("%s",s);  
	        find(s);  
	        cout<<num<<endl;  
	    }
		return 0;  
	} 
```



#### **树状数组维护fail树的dfs序**

给出n个字符串，表示n个人名，有两种操作：

?string ，统计字符串string中出现的属于城市居民的次数。

+id，把编号为id的人变为城市居民，如果已经是忽略。

-id，把编号为id的人变为不是城市居民，如果已经不是的话忽略。

现有m个操作，对于?输出结果。

```c++
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
#define lowbit(x) (x&-x)
#define maxn 1000010  
int n,m,tot;  
char s[maxn];  
vector<int> vec[maxn]; 								
int ch[maxn][26],fail[maxn],val[maxn],last[maxn]; 
int c[maxn],in[maxn],out[maxn],tim,id[maxn],use[maxn];
inline void Modify(int x,int num)
{
	if(x==0) return ;
	while(x<maxn)
		c[x]+=num,x+=lowbit(x);
}
inline void Add(int x,int y,int num)
{
	Modify(x,num);Modify(y,-num);	
}
inline int Query(int x)
{
	int res=0;
	while(x>0)
		res+=c[x],x-=lowbit(x);	
	return res;
}
inline void Init()  
{  
	tot=1;tim=0;
    memset(ch[0],0,sizeof(ch[0]));  
	memset(val,0,sizeof(val)); 
	memset(use,0,sizeof(use)); 
}     
inline int idx(char c){ return c-'a';}  
void Insert(char*s,int x)
{  
    int u=0,len=strlen(s);  
    for(int i=0;i<len;++i)  
    {  
        int c=idx(s[i]);  
        if(!ch[u][c])  
        {  
            memset(ch[tot],0,sizeof(ch[tot]));  
            val[tot]=0;  
            ch[u][c]=tot++;   
        }  
        u=ch[u][c];  
    }  
    val[u]=x;
    id[x]=u;
}
void GetFail()  
{  
    queue<int> q;  
    fail[0]=0;
    for(int c=0;c<26;++c)  
    {  
        int u=ch[0][c];  
        if(u){ fail[u]=0;q.push(u);last[u]=0; }  
    }  
    while(!q.empty())  
    {  
        int r=q.front(); q.pop();  
        vec[fail[r]].push_back(r);
        for(int c=0;c<26;++c)  
        {  								
            int u=ch[r][c];  
            if(!u){ch[r][c]=ch[fail[r]][c];continue;}  
            q.push(u);  
            int v=fail[r];  
            fail[u]=ch[v][c];  
            last[u] = val[fail[u]]?fail[u]:last[fail[u]];  
        }  
    }  
} 
void dfs(int u)
{
	in[u]=++tim;
	for(int i=0,len=vec[u].size();i<len;++i)
		dfs(vec[u][i]);
	out[u]=tim;
}
void clac(int x,int num)
{
	Add(in[id[x]],out[id[x]]+1,num);
}
void work()
{
	ll ans=0;
	int u=0,len=strlen(s);
	for(int i=1;i<len;++i)
	{
		int r=idx(s[i]);
		u=ch[u][r];
		ans+=Query(in[u]);	
	}
	printf("%lld\n",ans);
}
int main()
{
	scanf("%d%d",&m,&n);
	Init();
	for(int i=1;i<=n;++i) 
		scanf("%s",s),Insert(s,i),use[i]=1;
	GetFail();
	dfs(0);
	for(int i=1;i<=n;++i) clac(i,1);
 	while(m--)
	{
        scanf("%s",s);
        if(s[0]=='?') work();
        else
		{
			int x;
            sscanf(s+1,"%d",&x);
            if(use[x]&&s[0]=='-')
			{
                use[x] = 0;
                clac(x,-1);
            }
			else if(!use[x]&&s[0]=='+')
			{
                use[x] = 1;
                clac(x,1);
            }
        }
    }
	return 0;
}
```



#### **主席树维护fail树的dfs序**

给你n个单词，Q个询问: l ,r ,x在第l到r个单词之间，第x个单词出现了多少次。CF547E(CSL)

```c++
using namespace std;
typedef long long ll;
const int INF=0x3f3f3f3f;
const int maxn=2e5+10;
int N,Q;
char s[maxn];
int sz,ch[maxn][26],fail[maxn],val[maxn];
int in[maxn],out[maxn],tim,Len[maxn],tot,T[maxn];
vector<int> vec[maxn],ts[maxn];
void Init()
{
    sz=1; tim=tot=0;
    memset(ch[sz],0,sizeof ch[sz]);
    memset(val,0,sizeof val);
}
int idx(char c){return c-'a';}
void Insert(char *s,int x)
{
    Len[x]=strlen(s);int u=0;
    for(int i=0;i<Len[x];++i)
    {
        int c=idx(s[i]); ts[x].push_back(c);
        if(!ch[u][c]){memset(ch[sz],0,sizeof ch[sz]);ch[u][c]=sz++;}
        u=ch[u][c];
    }
    val[x]=u;
}
void GetFail()
{
    queue<int> q;
    fail[0]=0;
    for(int i=0;i<26;++i)
        if(ch[0][i]) {fail[ch[0][i]]=0;q.push(ch[0][i]);}
    while(!q.empty())
    {
        int u=q.front();q.pop();
        vec[fail[u]].push_back(u);
        for(int i=0;i<26;++i)
        {
            int c=ch[u][i];
            if(!c){ch[u][i]=ch[fail[u]][i];continue;}
            q.push(c);
            fail[c]=ch[fail[u]][i];
        }
    }
}
void dfs(int u)
{
    in[u]=++tim;
    for(int i=0,len=vec[u].size();i<len;++i)
        dfs(vec[u][i]);
    out[u]=tim;
}

struct Tree{
    int ls,rs;
    int num;
} tr[maxn*40];
inline void Update(int y,int &x,int l,int r,int pos)
{
    tr[++tot]=tr[y];tr[tot].num++; x=tot;
    if(l==r) return;
    int mid=l+r>>1;
    if(pos<=mid) Update(tr[y].ls,tr[x].ls,l,mid,pos);
    else Update(tr[y].rs,tr[x].rs,mid+1,r,pos);
}
inline void Build(int y,int &x,int id)
{
    int now=0,last=y;
    for(int i=0;i<Len[id];++i)
    {
        now=ch[now][ts[id][i]];
        Update(last,x,1,tim,in[now]);
        last=x;
    }
}
inline int Query(int y,int x,int l,int r,int L,int R)
{
    if(L<=l&&r<=R) return tr[x].num-tr[y].num;
    int mid=l+r>>1,res=0;
    if(L<=mid) res+=Query(tr[y].ls,tr[x].ls,l,mid,L,R);
    if(R>mid) res+=Query(tr[y].rs,tr[x].rs,mid+1,r,L,R);
    return res;
}

int main()
{
    scanf("%d%d",&N,&Q);
    Init();
    for(int i=1;i<=N;++i) scanf("%s",s),Insert(s,i);
    GetFail();
    dfs(0);
    for(int i=1;i<=N;++i) Build(T[i-1],T[i],i);
    while(Q--)
    {
        int l,r,x;
        scanf("%d%d%d",&l,&r,&x);
        printf("%d\n",Query(T[l-1],T[r],1,tim,in[val[x]],out[val[x]]));
    }

    return 0;
}
```



#### **长度不超过m的串所得到的最大权值**

fail树上状压DP

n个串，每个串一个值，问你在最长为m的字符串能得到的最大权值。(重复只算一次)

```c++
#include<bits/stdc++.h>
using namespace std;
const int INF=0x3f3f3f3f;
const int maxn=1010;
int ch[maxn][4],val[maxn],fail[maxn];
int w[12],tot,n,m,Ans[maxn];
bool dp[2][maxn][1<<11];
char s[110];

void Init()
{
	tot=1;
	memset(val,0,sizeof val);
	memset(ch[tot],0,sizeof ch[tot]);	
}

void Insert(char *s,int x)
{
	int len=strlen(s),u=0;
	for(int i=0;i<len;++i)
	{
		int c=s[i]-'a';
		if(!ch[u][c]) {ch[u][c]=tot++;memset(ch[tot],0,sizeof ch[tot]);}
		u=ch[u][c];
	}
	val[u]=1<<x;
}

void GetFail()
{
	queue<int>q;
	for(int i=0;i<4;++i)
		if(ch[0][i]) q.push(ch[0][i]);
	while(!q.empty())
	{
		int u=q.front();q.pop();
		for(int i=0;i<4;++i)
		{
			int v=ch[u][i];
			if(!v){ch[u][i]=ch[fail[u]][i];continue;}
			q.push(v);
			fail[v]=ch[fail[u]][i];
			val[v]|=val[fail[v]];
		}
	}
}

void Work()
{
	dp[0][0][0]=1;
	for(int i=1;i<=m;++i)
	{
		memset(dp[i&1],0,sizeof dp[i&1]);
		for(int j=0;j<tot;++j)
			for(int k=0;k<4;++k)
				for(int z=0;z<(1<<n);++z)
				{
					if(dp[(i-1)&1][j][z])
						dp[i&1][ch[j][k]][z|val[ch[j][k]]]=1;	
				}
	}
}

int GetAns(int x)
{
	int ans=0;
	for(int i=0;i<n;++i)
		if(x&(1<<i)) ans+=w[i];
	return ans;	
}

int main()
{
	scanf("%d%d",&n,&m);
	Init();
	for(int i=0;i<n;++i)
		scanf("%s%d",s,w+i),Insert(s,i);
	GetFail();
	Work();
	int res=-INF;
	
	for(int j=0;j<(1<<n);++j) Ans[j]=GetAns(j);
	for(int i=0;i<tot;++i)
		for(int j=0;j<(1<<n);++j)
			if(dp[m&1][i][j]) res=max(res,Ans[j]);
	
	if(res<0) puts("Unhappy!");
	else printf("%d\n",res);
	
	return 0;	
}
```



#### [**DP+AC自动机+最短路**](https://www.cnblogs.com/Konjakmoyu/p/5665802.html)

给出n个资源，m个病毒，将资源串拼接成一个串，必须包含所有的资源串，可以重叠，但是不能包含病毒。问最小的长度为多少。

```c++
const int maxn=1010;
const int maxnode=60010;
char s[maxn];
int n,m,tot,ans,ln[20];
int pre[maxn],vis[maxnode],dis[20][maxnode],dp[(1<<10)][15];
int ch[maxnode][2],link[maxnode],fail[maxnode],val[maxnode];
void Init()
{ 
    tot=1; ans=INF; 
    memset(ch[0],0,sizeof(ch[0]));
}
int GetId(char ch){return ch-'0';}
void Insert(int tk)
{
    scanf("%s",s);
    int len=strlen(s),u=0;
    if(tk!=-1) ln[tk]=len;
    for(int i=0;i<len;++i)
    {
        int x=GetId(s[i]);
        if(!ch[u][x])
        {
            memset(ch[tot],0,sizeof(ch[tot]));
            val[tot]=0;
            ch[u][x]=tot++;
        }
        u=ch[u][x];
    }
    if(tk==-1) val[u]=-1;
    else val[u]=tk,pre[tk]=u;
}
void GetFail()
{
    queue<int> q; 
    q.push(0); fail[0]=0;
    while(!q.empty())
    {
        int u=q.front();q.pop();
        for(int i=0;i<=1;++i)
        {
            if(ch[u][i])
            {
                fail[ch[u][i]] = u?ch[fail[u]][i]:0;
                q.push(ch[u][i]);
            }
            else ch[u][i]=ch[fail[u]][i];
        }
    }
}
void spfa(int s)
{
    memset(dis[s],INF,sizeof(dis[s]));
    memset(vis,0,sizeof(vis));
    dis[s][pre[s]]=0;
    queue<int> q;
    q.push(pre[s]);vis[pre[s]]=1;
    while(!q.empty())
    {
        int u=q.front();
        for(int i=0;i<=1;++i)
        {
            if(val[ch[u][i]]!=-1)
            {
                if(dis[s][ch[u][i]]>dis[s][u]+1)
                {
                    dis[s][ch[u][i]]=dis[s][u]+1;
                    if(!vis[ch[u][i]])
                    {
                        vis[ch[u][i]]=1;
                        q.push(ch[u][i]);
                    }
                }
            }
        } 
        q.pop(); vis[u]=0;
    }
}
bool check(int x,int y,int z)
{
    if(((1<<y-1)&z)==0) return false;
    int num=0;
    for(int i=1;i<=n;++i){ if((1<<i-1)&z) num++; }
    return num==x? true:false;
}
void work()
{
    memset(dp,INF,sizeof(dp));
    for(int i=1;i<=n;++i) dp[(1<<i-1)][i]=ln[i];
    for(int i=1;i<n;++i)
    {
        for(int j=1;j<=n;++j)
        {
            for(int k=1;k<=(1<<n)-1;++k)
            {
                if(check(i,j,k))
                {
                    for(int l=1;l<=n;++l)
                    {
                        if((k&(1<<l-1)) == 0) 
                            dp[k+(1<<l-1)][l]=min(dp[k+(1<<l-1)][l],dp[k][j]+dis[j][pre[l]]);
                    }
                }
            } 
        }
    }
    for(int i=1;i<=n;++i) ans=min(ans,dp[(1<<n)-1][i]);
}
int main()
{
    while(scanf("%d%d",&n,&m)&&n&&m)
    {
        Init();
        for(int i=1;i<=n;++i) Insert(i);
        for(int i=1;i<=m;++i) Insert(-1);
        GetFail();
        for(int i=1;i<=n;++i) spfa(i);
        work();
        printf("%d\n",ans);    
    }
    return 0;
}
```



### **PAM(回文自动机)**

#### **模板**

```c++
struct Palindromic_Tree{
	int next[MAXN][26];//next指针，next指针和字典树类似，指向的串为当前串两端加上同一个字符构成 
	int fail[MAXN];;//fail指针，失配后跳转到fail指针指向的节点
	int cnt[MAXN];//表示i表示的回文字符串在整个字符串中出现了多少次（建树时求出的不是完全的，最后count()函数跑一遍以后才是正确的）
	int num[MAXN];//表示以节点i表示的最长回文串的最右端点为回文串结尾的回文串个数
	int len[MAXN];//len[i]表示节点i表示的回文串的长度
	int S[MAXN];//存放添加的字符
	int last;//指向上一个字符所在的节点，方便下一次add
	int n;//字符数组指针表示添加的字符个数。
	int p;//节点指针,表示添加的节点个数。(p-2):表示本质不同的回文串数量 
 
	int newnode(int l) 
	{
		for(int i=0;i<26;++i) next[p][i]=0;
		cnt[p]=0;
		num[p]=0;
		len[p]=l;
		return p++;
	}
 
	void Init() 
	{
		p=0;
		newnode( 0);
		newnode(-1);
		last=0;
		n=0;
		S[n]=-1;
		fail[0]=1;
	}
 
	int get_fail(int x)
	{
		while(S[n-len[x]-1]!=S[n])x=fail[x] ;
		return x ;
	}
 
	void add(int c) 
	{
		S[++ n]=c;
		int cur=get_fail(last) ;
		if(!next[cur][c]) 
		{
			int now=newnode(len[cur]+2) ;
			fail[now]=next[get_fail(fail[cur])][c] ;
			next[cur][c]=now ;
			num[now]=num[fail[now]]+1;
		}
		last=next[cur][c];
		cnt[last]++;
	}
 
	ll count()  
	{
		ll res=p*1ll;
		for(int i=p-1;i>=0;--i) cnt[fail[i]]+=cnt[i];
//父亲累加儿子的cnt，因为如果fail[v]=u，则u一定是v的子回文串！
		return (res-2);//本质不同的回文串数量 
	}
} pam;

```



#### **一个串里所有本质不同的回文子串满足一个串是另一个的子串****的对数**

```c++
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
const int maxn=1e5+5;
int T;
char s[maxn];
ll ans;
struct PAM{
    int nxt[maxn][26],fail[maxn],len[maxn];
    int S[maxn],dp[maxn],vis[maxn];
    int last,n,p;
    int NewNode(int l) 
	{
        memset(nxt[p],0,sizeof(nxt[p]));
        len[p]=l;
        dp[p]=0;
        return p++;
    }
    void Init() 
	{
        ans=0; n=last=p=0;
        NewNode(0); NewNode(-1);
        S[n]=-1; fail[0]=1;
    }
    int GetFail(int x) 
	{
        while(S[n-len[x]-1]!=S[n]) x=fail[x];
        return x;
    }
    void add(int c) 
	{
        S[++n]=c;
        int cur=GetFail(last);
        if(!nxt[cur][c]) 
		{
            int now=NewNode(len[cur]+2);
            fail[now]=nxt[GetFail(fail[cur])][c];
            nxt[cur][c]=now;
        }
        last=nxt[cur][c];
    }
    void dfs(int x,int fa) 
	{
		int cnt=0,cx=x; vis[x]=1;
        while(fail[cx]!=0&&fail[cx]!=1&&!vis[fail[cx]]) 
            cx=fail[cx],vis[cx]=1,++cnt;
        dp[x]=cnt;
        if(x!=1&&x!=0&&fa!=0&&fa!=1) 
            dp[x]=dp[fa]+cnt+1;
        ans+=dp[x];
        for(int i=0;i<26;++i){if(nxt[x][i]) dfs(nxt[x][i],x);}
        vis[x]=0; cx=x;
        while(cnt--) cx=fail[cx],vis[cx]=0;
    }
    void work()
    {
    	Init();
        for(int i=1,len=strlen(s+1);i<=len;i++) add(s[i]-'a');
    	dfs(1,1);
    	dfs(0,0);
    }

} pam;
int main() 
{
    scanf("%d", &T);
    for(int cas=1;cas<=T;++cas) 
	{
        scanf("%s",s+1);
        pam.work();
        printf("Case #%d: %lld\n",cas,ans);
    }
    return 0;
}
```

 

#### **求公共回文串个数** 

```c++
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
const int maxn=2e5+100;
struct PAM{
	int nxt[maxn][26],len[maxn],cnt[maxn],fail[maxn];
	int S[maxn],last,p,now;
	int NewNode(int l)
	{
		memset(nxt[p],0,sizeof nxt[p]);
		cnt[p]=0;len[p]=l;
		return p++;
	}
	void Init()
	{
		p=0;
		NewNode(0);NewNode(-1);
		last=0;now=0;
		S[now++]=-1;fail[0]=1;
	}
	int GetFail(int x)
	{
		while(S[now-len[x]-2]!=S[now-1]) x=fail[x];
		return x;
	} 
	void Add(int c)
	{
		S[now++]=c;
		int cur=GetFail(last);
		if(!nxt[cur][c])
		{
			int now=NewNode(len[cur]+2);
			fail[now]=nxt[GetFail(fail[cur])][c];
			nxt[cur][c]=now;
		}
		last=nxt[cur][c];cnt[last]++;
	}
	void Count()
	{
		for(int i=p-1;i>=0;i--) cnt[fail[i]]+=cnt[i];
		cnt[0]=cnt[1]=0;
	}
 
} pam1,pam2;

ll dfs(int u,int v)
{
	ll res=0;
	for (int i=0;i<26;i++)
	{
		int uu=pam1.nxt[u][i];
		int vv=pam2.nxt[v][i];
		if(uu&&vv)
		{
			res+=1LL*pam1.cnt[uu]*pam2.cnt[vv];
			res+=dfs(uu,vv);
		}
	}
	return res;
}
int T,len1,len2;
char s1[maxn],s2[maxn];
int main()
{
	scanf("%d",&T);
	for(int cas=1;cas<=T;++cas)
	{
		pam1.Init(); pam2.Init();
		scanf("%s%s",s1,s2);
		len1=strlen(s1); len2=strlen(s2);
		for(int i=0;i<len1;i++) pam1.Add(s1[i]-'a');
		for(int i=0;i<len2;i++) pam2.Add(s2[i]-'a');
		pam1.Count();  pam2.Count();
		printf("Case #%d: %lld\n",cas++,dfs(0,0)+dfs(1,1));
	}
	return 0;
}
```

 

#### **邻接表优化求 相交回文串对 个数**

```c++
const int N=2e6+5;
const int mod=51123987;
int n,fa[N],len[N],dep[N],tot,last,p1[N],p2[N],ans;
int to[N],nxt[N],ww[N],head[N],cnt;
char s[N];
void init()
{
    fa[last=0]=fa[1]=1;
    len[0]=0;len[tot=1]=-1;
    memset(head,0,sizeof(head));cnt=0;
}
void link(int u,int v,int c)
{
    to[++cnt]=v;nxt[cnt]=head[u];ww[cnt]=c;
    head[u]=cnt;
}
int tr(int v,int c)
{
    for(int e=head[v];e;e=nxt[e])
        if(ww[e]==c) return to[e];
    return 0;
}
void extend(int c,int n)
{
    int v=last;
    while(s[n-len[v]-1]!=s[n]) v=fa[v];
    if(!tr(v,c))
    {
        int u=++tot,k=fa[v];
        len[u]=len[v]+2;
        while(s[n-len[k]-1]!=s[n]) k=fa[k];
        fa[u]=tr(k,c); dep[u]=dep[fa[u]]+1;
        link(v,u,c);
    }
    last=tr(v,c);
}
int main()
{
    scanf("%d",&n);
    scanf("%s",s+1);
    init();
    for(int i=1;i<=n;++i) extend(s[i]-'a',i),(ans+=(p1[i]=dep[last]))%=mod;
    ans=1ll*ans*(ans-1)/2%mod;
    reverse(s+1,s+n+1);
    init();
    for(int i=1;i<=n;++i) extend(s[i]-'a',i),p2[n-i+1]=dep[last];
    for(int i=n;i;--i) (p2[i]+=p2[i+1])%=mod;
    for(int i=1;i<=n;++i) ans=(ans-1ll*p1[i]*p2[i+1]%mod+mod)%mod;
    printf("%d\n",ans);
    return 0;
}

```



### **序列自动机**

next[i][j]:表示在原串s第i位后面的第一个j出现的位置 

```c++
for(int i=n;i;i--)
{
	for(int j=1;j<=a;j++) next[i-1][j]=next[i][j];
	next[i-1][s[i]]=i;
}
```



#### **求子序列个数**

```c++
int dfs(int x)//0
{
	if(f[x]) return f[x];
	for(r int i=1;i<=a;i++)
		if(next[x][i]) f[x]+=dfs(next[x][i]);
	return ++f[x];
}
```



#### **求两个串的公共子序列个数**

```c++
int dfs(int x,int y)//0 0
{//表示目前字符是串1的第x位,串2的第y位
	if(f[x][y]) return f[x][y];
	for(r int i=1;i<=a;i++)
		if(next1[x][i]&&next2[y][i]) f[x][y]+=dfs(next1[x][i],next2[y][i]);
	return ++f[x][y];
}
```



#### **求串回文子序列个数**

```c++
int dfs(int x,int y)//算空串
{
	if(f[x][y]) return f[x][y];
	for(int i=1;i<=26;i++)
		if(nxt1[x][i]&&nxt2[y][i])
		{
			if(nxt1[x][i]+nxt2[y][i]>n+1) continue;//x+y>n+1,状态不合法,不进行统计
			if(nxt1[x][i]+nxt2[y][i]<n+1) f[x][y]++;
			//满足x+y=n+1的奇串不会被漏掉,其他奇串需要特别统计
			f[x][y]=(f[x][y]+dfs(nxt1[x][i],nxt2[y][i]));
		}
	return ++f[x][y];
}
```



#### **求A,B的最长公共子序列S使得C是S的子序列**

还是同样的Dfs(LLx,LLy,LLz)，表示一匹配到C的z位
改变一下C的构建方法

```c++
for(LL i=1;i<=a;++i) nxt[n][i]=n;
for(LL i=0;i<n;++i)
{
    for(LL j=1;j<=a;++j) nxt[i][j]=i;
    nxt[i][c[i+1]]=i+1;
}

```

>  