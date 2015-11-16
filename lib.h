#define buf_length 100

void init();
double convertStringToDouble(char *str);
void creatVoltList(FILE *fp);
void creatAmberList(FILE *fp);
void creatResistanceList(FILE *fp);
void creatCapacitorList(FILE *fp);
void creatInductorList(FILE *fp);
void creatDiodeList(FILE *fp);
void creatMOSList(FILE *fp);
void creatBJTList(FILE *fp);
void analysis(FILE *fp);

typedef struct Volt{
  char *name;
  int node1;
  int node2;
  double value;
  struct Volt *next;
}VoltT;

typedef struct Amper{
  char *name;
  int node1;
  int node2;
  double value;
  struct Amper *next;
}AmperT;

typedef struct Resistance{
  char *name;
  int node1;
  int node2;
  double value;
  struct Resistance *next;  
}ResistanceT;

typedef struct Capacitor{
  char *name;
  int node1;
  int node2;
  double value;
  struct Capacitor *next;
}CapacitorT;

typedef struct Inductor{
  char *name;
  int node1;
  int node2;
  double value;
  struct Inductor *next;
}InductorT;

typedef struct Diode{
  char *name;
  int node1;
  int node2;
  struct Diode *next;
}DiodeT;

typedef struct Mos{
  char *name;
  int D;
  int G;
  int S;
  struct Mos *next;
}MosT;

typedef struct Bjt{
  char *name;
  int C;
  int B;
  int E;
  struct Bjt *next;  
}BjtT;

typedef struct nodes{
 int name;
 struct nodes *next;
}node;

//Dilwsi twn katholikwn metavlitwn tou programmatos (roots kai elegxou gia geiwsi)
VoltT *rootV;
AmperT *rootI;
ResistanceT *rootR;
CapacitorT *rootC;
InductorT *rootL;
DiodeT *rootD;
MosT *rootM;
BjtT *rootB;
int ground;
node *root_node, *pre_node;
int use_lu, use_cholesky,found_dc_sweep,source,found_iter, sparse_option, sweep_node1, sweep_node2, sweep_value;
float start_value,end_value,step,itol_value;

VoltT *currV;
AmperT *currI;
InductorT *currL;
ResistanceT *currR;
node *curr_node;

int count_nodes();
void printListOf2Node(void* list);
void printNodeList();
void printVoltList(VoltT *list);
void printAmperList(AmperT *list);
void printResistanceList(ResistanceT *list);
void printCapacitorList(CapacitorT *list);
void printInductorList(InductorT *list);
void printDiodeList(DiodeT *list);
void printMosList(MosT *list);
void printBjttList(BjtT *list);

