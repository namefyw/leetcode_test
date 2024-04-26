
#include "leetcode.h"

#if leetcode_20
char pairs(char a) {
    if (a == '}') return '{';
    if (a == ']') return '[';
    if (a == ')') return '(';
    return 0;
}

bool isValid(char* s) {
    int n = strlen(s);
    if (n % 2 == 1) {
        return false;
    }
    int stk[n + 1], top = 0;
    for (int i = 0; i < n; i++) {
        char ch = pairs(s[i]);
        if (ch) {
            if (top == 0 || stk[top - 1] != ch) {
                return false;
            }
            top--;
        } else {
            stk[top++] = s[i];
        }
    }
    return top == 0;
}

#endif

#if leetcode_21
/**
 * Definition for singly-linked list.
 * struct ListNode {
 *     int val;
 *     struct ListNode *next;
 * };
 */
#if 1 //实体指针法
struct ListNode* mergeTwoLists(struct ListNode* list1, struct ListNode* list2) {
    if (list1 == NULL) {
        return list2;
    }
    if (list2 == NULL) {
        return list1;
    }
    struct ListNode head;
    struct ListNode* tmp = &head;
    while (list1 && list2) {
        if (list1->val < list2->val) {
            tmp->next = list1;
            list1 = list1->next;
        } else {
            tmp->next = list2;
            list2 = list2->next;
        }
        tmp = tmp->next;
    }
    tmp->next = (list1 != NULL) ? list1 : list2;
    return head.next;
}
#endif
#endif

#if leetcode_219
typedef struct record{
  int data;
  int index;
}record;

int cmp(const void*a, const void*b){
  return (*(record*)a).data - (*(record*)b).data;
}
bool containsNearbyDuplicate(int* nums, int numsSize, int k){
  if(numsSize == 0 || numsSize == 1){
      return false;
  }
  record temp[numsSize];
  for(int i = 0; i < numsSize; i++){
      temp[i].data = nums[i];
      temp[i].index = i;
  }
  qsort(temp,numsSize,sizeof(record),cmp);
  for(int i = 0; i < numsSize - 1; i++){
      if(temp[i].data == temp[i+1].data){
          if(abs(temp[i].index - temp[i+1].index) <= k){
              return true;
          }
      }
  }
  return false;
}
#endif

#if leetcode_1793
int maximumScore(int* nums, int numsSize, int k){
    int ans = nums[k];
    int left = k;
    int right = k;
    int min = nums[k];
    while (left > 0 || right < numsSize - 1) {
        if (right >= numsSize - 1 || (left > 0 && nums[left - 1] > nums[right + 1])) {
            left--;
            min = fmin(min, nums[left]);
        } else {
            right++;
            min = fmin(min, nums[right]);
        }
        ans = fmax(ans, (right - left + 1) * min);
    }
    return ans;
}
#endif

#if leetcode_2312
#if 1
typedef struct {
    long long key;
    long long val;
    UT_hash_handle hh;
} HashItem;

HashItem *hashFindItem(HashItem **obj, long long key) {
    HashItem *pEntry = NULL;
    HASH_FIND(hh, *obj, &key, sizeof(long long), pEntry);
    return pEntry;
}

bool hashAddItem(HashItem **obj, long long key, long long val) {
    if (hashFindItem(obj, key)) {
        return false;
    }
    HashItem *pEntry = (HashItem *)malloc(sizeof(HashItem));
    pEntry->key = key;
    pEntry->val = val;
    HASH_ADD(hh, *obj, key, sizeof(key), pEntry);
    return true;
}

bool hashSetItem(HashItem **obj, long long key, long long val) {
    HashItem *pEntry = hashFindItem(obj, key);
    if (!pEntry) {
        hashAddItem(obj, key, val);
    } else {
        pEntry->val = val;
    }
    return true;
}

int hashGetItem(HashItem **obj, long long key, long long defaultVal) {
    HashItem *pEntry = hashFindItem(obj, key);
    if (!pEntry) {
        return defaultVal;
    }
    return pEntry->val;
}

void hashFree(HashItem **obj) {
    HashItem *curr = NULL, *tmp = NULL;
    HASH_ITER(hh, *obj, curr, tmp) {
        HASH_DEL(*obj, curr);
        free(curr);
    }
}

long long pairHash(int x, int y) {
    return ((long long) x << 16) ^ y;
}

long long dfs(int x, int y, long long **memo, HashItem **val) {
    if (memo[x][y] != -1) {
        return memo[x][y];
    }

    long long ret = hashGetItem(val, pairHash(x, y), 0LL);
    if (x > 1) {
        for (int i = 1; i < x; ++i) {
            ret = fmax(ret, dfs(i, y, memo, val) + dfs(x - i, y, memo, val));
        }
    }
    if (y > 1) {
        for (int j = 1; j < y; ++j) {
            ret = fmax(ret, dfs(x, j, memo, val) + dfs(x, y - j, memo, val));
        }
    }
    return memo[x][y] = ret;
};

long long sellingWood(int m, int n, int** prices, int pricesSize, int* pricesColSize) {
    HashItem *value = NULL;
    long long *memo[m + 1];
    for (int i = 0; i <= m; i++) {
        memo[i] = (long long *)malloc(sizeof(long long) * (n + 1));
        for (int j = 0; j <= n; j++) {
            memo[i][j] = -1;
        }
    }
    for (int i = 0; i < pricesSize; ++i) {
        hashAddItem(&value, pairHash(prices[i][0], prices[i][1]), prices[i][2]);
    }
    long long ret = dfs(m, n, memo, &value);
    for (int i = 0; i <= m; i++) {
        free(memo[i]);
    }
    hashFree(&value);
    return ret;
}
#else
long long sellingWood(int m, int n, int** prices, int pricesSize, int* pricesColSize)
{
    long long int dp[m + 1][n + 1];
    int i, j;
    int col,row;
    long long int max,tmp;

    memset(dp, 0, sizeof(dp));

    // 初始化，将个块分割填入
    for (i = 0; i < pricesSize; ++i) {
        dp[prices[i][0]][prices[i][1]] = dp[prices[i][0]][prices[i][1]] > prices[i][2] ?
            dp[prices[i][0]][prices[i][1]] : prices[i][2];
    }
    // 遍历二维数组
    for (i = 1; i <= m ; ++i) {
        for (j = 1; j <= n; j++) {
            max = 0;
            // 横向分割;i/2，优化重复计算
            for (col = 1; col <= i / 2; col++) {
                tmp = dp[col][j] + dp[i-col][j];
                max = max > tmp? max : tmp;
            }
            // 纵向分割
            for (row = 1; row <= j/ 2; row++){
                tmp = dp[i][row] + dp[i][j-row];
                max = max > tmp? max : tmp;
            }
            dp[i][j] = dp[i][j] > max? dp[i][j] : max;
        }
    }
    return dp[m][n];
}
#endif
#endif
#if (leetcode_2007)
#if 0
typedef struct {
    int key;
    int val;
    UT_hash_handle hh;
} HashItem;

HashItem *hashFindItem(HashItem **obj, int key) {
    HashItem *pEntry = NULL;
    HASH_FIND_INT(*obj, &key, pEntry);
    return pEntry;
}

bool hashAddItem(HashItem **obj, int key, int val) {
    if (hashFindItem(obj, key)) {
        return false;
    }
    HashItem *pEntry = (HashItem *)malloc(sizeof(HashItem));
    pEntry->key = key;
    pEntry->val = val;
    HASH_ADD_INT(*obj, key, pEntry);
    return true;
}

bool hashSetItem(HashItem **obj, int key, int val) {
    HashItem *pEntry = hashFindItem(obj, key);
    if (!pEntry) {
        hashAddItem(obj, key, val);
    } else {
        pEntry->val = val;
    }
    return true;
}

int hashGetItem(HashItem **obj, int key, int defaultVal) {
    HashItem *pEntry = hashFindItem(obj, key);
    if (!pEntry) {
        return defaultVal;
    }
    return pEntry->val;
}

void hashFree(HashItem **obj) {
    HashItem *curr = NULL, *tmp = NULL;
    HASH_ITER(hh, *obj, curr, tmp) {
        HASH_DEL(*obj, curr);
        free(curr);
    }
}

static int cmp(const void *a, const void *b) {
    return *(int *)a - *(int *)b;
}

int* findOriginalArray(int* changed, int changedSize, int* returnSize) {
    qsort(changed, changedSize, sizeof(int), cmp);
    HashItem *count = NULL;
    for (int i = 0; i < changedSize; i++) {
        int a = changed[i];
        hashSetItem(&count, a, hashGetItem(&count, a, 0) + 1);
    }

    int *res = (int *)malloc(sizeof(int) * changedSize);
    int pos = 0;
    for (int i = 0; i < changedSize; i++) {
        int a = changed[i];
        if (hashGetItem(&count, a, 0) == 0) {
            continue;
        }
        hashSetItem(&count, a, hashGetItem(&count, a, 0) - 1);
        if (hashGetItem(&count, 2 * a, 0) == 0) {
            hashFree(&count);
            *returnSize = 0;
            return NULL;
        }
        hashSetItem(&count, a * 2, hashGetItem(&count, a * 2, 0) - 1);
        res[pos++] = a;
    }
    *returnSize = pos;
    hashFree(&count);
    return res;
}
#else
int cmp(const void* a, const void* b) {
    return *(int*) a - *(int*) b;
}

int* findOriginalArray(int* changed, int changedSize, int* returnSize) {
    qsort(changed, changedSize, sizeof(int), cmp); //排序的函数
    int* ans = malloc(changedSize / 2 * sizeof(int));
    int ans_idx = 0;
    int* q = malloc(changedSize / 2 * sizeof(int));
    int front = 0, rear = 0;
    for (int i = 0; i < changedSize; ++i) {
        int x = changed[i];
        if (front < rear) {
            if (q[front] < x) { // 无法配对
                free(ans);
                free(q);
                *returnSize = 0;
                return NULL;
            }
            if (q[front] == x) { // 配对成功
                ++front; // 清除一个标记
                continue;
            }
        }
        if (ans_idx == changedSize / 2) {
            free(ans);
            free(q);
            *returnSize = 0;
            return NULL;
        }
        ans[ans_idx++] = x;
        q[rear++] = x * 2; // 添加双倍标记
    }
    *returnSize = ans_idx;
    return ans;
}
#endif
#endif
#if (leetcode_2642)
#if 1
/* 官方题解 */
typedef struct Edge {
    int to;
    int cost;
    struct Edge *next;
} Edge;

typedef struct {
    int first;
    int second;
} Node;

typedef bool (*cmp)(const void *, const void *);

typedef struct {
    Node *arr;
    int capacity;
    int queueSize;
    cmp compare;
} PriorityQueue;

Edge *createEdge(int to, int cost) {
    Edge *obj = (Edge *)malloc(sizeof(Edge));
    obj->to = to;
    obj->cost = cost;
    obj->next = NULL;
    return obj;
}

void freeEdgeList(Edge *list) {
    while (list) {
        Edge *p = list;
        list = list->next;
        free(p);
    }
}

Node *createNode(long long x, int y) {
    Node *obj = (Node *)malloc(sizeof(Node));
    obj->first = x;
    obj->second = y;
    return obj;
}

PriorityQueue *createPriorityQueue(int size, cmp compare) {
    PriorityQueue *obj = (PriorityQueue *)malloc(sizeof(PriorityQueue));
    obj->arr = (Node *)malloc(sizeof(Node) * size);
    obj->queueSize = 0;
    obj->capacity = size;
    obj->compare = compare;
    return obj;
}

static void swap(Node *arr, int i, int j) {
    Node tmp;
    memcpy(&tmp, &arr[i], sizeof(Node));
    memcpy(&arr[i], &arr[j], sizeof(Node));
    memcpy(&arr[j], &tmp, sizeof(Node));
}

static void down(Node *arr, int size, int i, cmp compare) {
    for (int k = 2 * i + 1; k < size; k = 2 * k + 1) {
        // 父节点 (k - 1) / 2，左子节点 k，右子节点 k + 1
        if (k + 1 < size && compare(&arr[k], &arr[k + 1])) {
            k++;
        }
        if (compare(&arr[k], &arr[(k - 1) / 2])) {
            break;
        }
        swap(arr, k, (k - 1) / 2);
    }
}

void Heapify(PriorityQueue *obj) {
    for (int i = obj->queueSize / 2 - 1; i >= 0; i--) {
        down(obj->arr, obj->queueSize, i, obj->compare);
    }
}

void Push(PriorityQueue *obj, Node *node) {
    memcpy(&obj->arr[obj->queueSize], node, sizeof(Node));
    for (int i = obj->queueSize; i > 0 && obj->compare(&obj->arr[(i - 1) / 2], &obj->arr[i]); i = (i - 1) / 2) {
        swap(obj->arr, i, (i - 1) / 2);
    }
    obj->queueSize++;
}

Node* Pop(PriorityQueue *obj) {
    swap(obj->arr, 0, obj->queueSize - 1);
    down(obj->arr, obj->queueSize - 1, 0, obj->compare);
    Node *ret =  &obj->arr[obj->queueSize - 1];
    obj->queueSize--;
    return ret;
}

bool isEmpty(PriorityQueue *obj) {
    return obj->queueSize == 0;
}

Node* Top(PriorityQueue *obj) {
    if (obj->queueSize == 0) {
        return NULL;
    } else {
        return &obj->arr[0];
    }
}

void FreePriorityQueue(PriorityQueue *obj) {
    free(obj->arr);
    free(obj);
}

bool greater(const void *a, const void *b) {
   return ((Node *)a)->first > ((Node *)b)->first;
}

typedef struct {
    Edge **graph;
    int nodeSize;
} Graph;


Graph* graphCreate(int n, int** edges, int edgesSize, int* edgesColSize) {
    Graph *obj = (Graph *)malloc(sizeof(Graph));
    obj->nodeSize = n;
    obj->graph = (Edge **)malloc(sizeof(Edge *) * n);
    for (int i = 0; i < n; i++) {
        obj->graph[i] = NULL;
    }
    for (int i = 0; i < edgesSize; i++) {
        int x = edges[i][0];
        int y = edges[i][1];
        int cost = edges[i][2];
        Edge *e = (Edge *)malloc(sizeof(Edge));
        e->to = y;
        e->cost = cost;
        e->next = obj->graph[x];
        obj->graph[x] = e;
    }

    return obj;
}

void graphAddEdge(Graph* obj, int* edge, int edgeSize) {
    int x = edge[0];
    int y = edge[1];
    int cost = edge[2];
    Edge *e = (Edge *)malloc(sizeof(Edge));
    e->to = y;
    e->cost = cost;
    e->next = obj->graph[x];
    obj->graph[x] = e;
}

int graphShortestPath(Graph* obj, int node1, int node2) {
    int n = obj->nodeSize;
    int dist[n];
    PriorityQueue *pq = createPriorityQueue(n * n, greater);
    for (int i = 0; i < n; i++) {
        dist[i] = INT_MAX;
    }
    dist[node1] = 0;
    Node val;
    val.first = 0;
    val.second = node1;
    Push(pq, &val);
    while (!isEmpty(pq)) {
        Node *p = Pop(pq);
        int cost = p->first;
        int cur = p->second;
        if (cur == node2) {
            return cost;
        }
        for (Edge *pEntry = obj->graph[cur]; pEntry; pEntry = pEntry->next) {
            int next = pEntry->to;
            int ncost = pEntry->cost;
            if (dist[next] > cost + ncost) {
                dist[next] = cost + ncost;
                val.first = cost + ncost;
                val.second = next;
                Push(pq, &val);
            }
        }
    }
    return -1;
}

void graphFree(Graph* obj) {
    for (int i = 0; i < obj->nodeSize; i++) {
        freeEdgeList(obj->graph[i]);
    }
    free(obj->graph);
    free(obj);
}
#endif
#if (q_err)
typedef struct {
    int vertices; // 顶点数量
    int **adjacency; // 邻接矩阵
} Graph;

// 创建图的实例
Graph* graphCreate(int n, int** edges, int edgesSize, int* edgesColSize) {
    Graph *graph = (Graph *)malloc(sizeof(Graph));
    graph->vertices = n;

    // 分配邻接矩阵的内存空间
    graph->adjacency = (int **)malloc(n * sizeof(int *));
    for (int i = 0; i < n; i++) {
        graph->adjacency[i] = (int *)calloc(n, sizeof(int));
    }

    // 根据边的信息构建邻接矩阵
    for (int i = 0; i < edgesSize; i++) {
        int src = edges[i][0];
        int dest = edges[i][1];
        int weight = edges[i][2];
        graph->adjacency[src][dest] = weight;
    }

    return graph;
}

// 向图中添加一条边
void graphAddEdge(Graph* obj, int* edge, int edgeSize) {
    int src = edge[0];
    int dest = edge[1];
    int weight = edge[2];
    obj->adjacency[src][dest] = weight;
}

// Dijkstra算法求最短路径
int graphShortestPath(Graph* obj, int node1, int node2) {
    int vertices = obj->vertices;
    int *dist = (int *)malloc(vertices * sizeof(int));
    int *visited = (int *)calloc(vertices, sizeof(int));

    // 初始化距离数组
    for (int i = 0; i < vertices; i++) {
        dist[i] = INT_MAX;
    }
    dist[node1] = 0;

    // 找到最短路径
    for (int count = 0; count < vertices - 1; count++) {
        int min_dist = INT_MAX;
        int min_index;

        // 选择未被访问的距离最小的顶点
        for (int v = 0; v < vertices; v++) {
            if (!visited[v] && dist[v] < min_dist) {
                min_dist = dist[v];
                min_index = v;
            }
        }

        visited[min_index] = 1; // 将选定的顶点标记为已访问

        // 更新距离数组
        for (int v = 0; v < vertices; v++) {
            if (!visited[v] && obj->adjacency[min_index][v] &&
                dist[min_index] != INT_MAX &&
                dist[min_index] + obj->adjacency[min_index][v] < dist[v]) {
                dist[v] = dist[min_index] + obj->adjacency[min_index][v];
            }
        }
    }

    int shortest_path = dist[node2];
    free(dist);
    free(visited);
    return shortest_path;
}

// 释放图的实例及其占用的内存
void graphFree(Graph* obj) {
    for (int i = 0; i < obj->vertices; i++) {
        free(obj->adjacency[i]);
    }
    free(obj->adjacency);
    free(obj);
}
#endif
#endif
#if (leetcode_2671)
#if 0
typedef struct {
    int key;
    int val;
    UT_hash_handle hh;
} HashItem;

HashItem *hashFindItem(HashItem **obj, int key) {
    HashItem *pEntry = NULL;
    HASH_FIND_INT(*obj, &key, pEntry);
    return pEntry;
}

bool hashAddItem(HashItem **obj, int key, int val) {
    if (hashFindItem(obj, key)) {
        return false;
    }
    HashItem *pEntry = (HashItem *)malloc(sizeof(HashItem));
    pEntry->key = key;
    pEntry->val = val;
    HASH_ADD_INT(*obj, key, pEntry);
    return true;
}

bool hashSetItem(HashItem **obj, int key, int val) {
    HashItem *pEntry = hashFindItem(obj, key);
    if (!pEntry) {
        hashAddItem(obj, key, val);
    } else {
        pEntry->val = val;
    }
    return true;
}

int hashGetItem(HashItem **obj, int key, int defaultVal) {
    HashItem *pEntry = hashFindItem(obj, key);
    if (!pEntry) {
        return defaultVal;
    }
    return pEntry->val;
}

void hashFree(HashItem **obj) {
    HashItem *curr = NULL, *tmp = NULL;
    HASH_ITER(hh, *obj, curr, tmp) {
        HASH_DEL(*obj, curr);
        free(curr);
    }
}

typedef struct {
    HashItem *freq;
    HashItem *freq_cnt;
} FrequencyTracker;


FrequencyTracker* frequencyTrackerCreate() {
    FrequencyTracker *obj = (FrequencyTracker *)malloc(sizeof(FrequencyTracker));
    obj->freq = NULL;
    obj->freq_cnt = NULL;
    return obj;
}

void frequencyTrackerAdd(FrequencyTracker* obj, int number) {
    int prev = hashGetItem(&obj->freq, number, 0);
    hashSetItem(&obj->freq_cnt, prev, hashGetItem(&obj->freq_cnt, prev, 0) - 1);
    hashSetItem(&obj->freq, number, prev + 1);
    hashSetItem(&obj->freq_cnt, prev + 1, hashGetItem(&obj->freq_cnt, prev + 1, 0) + 1);
}

void frequencyTrackerDeleteOne(FrequencyTracker* obj, int number) {
    int prev = hashGetItem(&obj->freq, number, 0);
    if (prev == 0) {
        return;
    }
    hashSetItem(&obj->freq_cnt, prev, hashGetItem(&obj->freq_cnt, prev, 0) - 1);
    hashSetItem(&obj->freq, number, prev - 1);
    hashSetItem(&obj->freq_cnt, prev - 1, hashGetItem(&obj->freq_cnt, prev - 1, 0) + 1);
}

bool frequencyTrackerHasFrequency(FrequencyTracker* obj, int frequency) {
    return hashGetItem(&obj->freq_cnt, frequency, 0) > 0;
}

void frequencyTrackerFree(FrequencyTracker* obj) {
    hashFree(&obj->freq);
    hashFree(&obj->freq_cnt);
    free(obj);
}
#endif
#if 1
#define N 100010
typedef struct {
    int *arr;
    int *freArr;
} FrequencyTracker;


FrequencyTracker* frequencyTrackerCreate() {
    FrequencyTracker* obj = (FrequencyTracker*)calloc(1, sizeof(FrequencyTracker));
    obj->arr = (int *)calloc(N, sizeof(int));
    obj->freArr = (int *)calloc(N, sizeof(int));
    return obj;
}

void frequencyTrackerAdd(FrequencyTracker* obj, int number) {
    obj->freArr[obj->arr[number]]--;
    obj->arr[number]++;
    obj->freArr[obj->arr[number]]++;
}

void frequencyTrackerDeleteOne(FrequencyTracker* obj, int number) {
    if(obj->arr[number] == 0) return;
    obj->freArr[obj->arr[number]]--;
    obj->arr[number]--;
    obj->freArr[obj->arr[number]]++;
}

bool frequencyTrackerHasFrequency(FrequencyTracker* obj, int frequency) {
    return obj->freArr[frequency];
}

void frequencyTrackerFree(FrequencyTracker* obj) {
    free(obj->arr);
    free(obj->freArr);
    free(obj);
}
#endif
#endif

