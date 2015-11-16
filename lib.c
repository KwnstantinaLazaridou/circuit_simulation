#include <stdio.h>
#include "lib.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lib_matrix_sparse.h"


void init(){
	rootV=NULL;
	rootI=NULL;
	rootR=NULL;
	rootC=NULL;
	rootL=NULL;
	rootD=NULL;
	rootM=NULL;
	rootB=NULL;
	ground=0;
        use_lu=1;
	use_cholesky=0;
        found_dc_sweep=0;
	found_iter=0;
	sparse_option=0;
	sparse_elements=0;
}

double convertStringToDouble(char *str){
    int i=0;
    int j=0;    
    char* str1=(char*)malloc(sizeof(char)*(i+1));
    double arith1=0;
    int arith2=0;

    for(i=0;i<strlen(str);i++)
    {
        if(str[i]=='e')
            break;
    }

    for(j=0;j<i;j++)
    {
        str1[j]=str[j];
    }
    str1[j]='\0';
    arith1=atof(str1);  
    free(str1);

    str1=(char*)malloc(sizeof(char)*(strlen(str)-i));
    i++;
    for(;i<strlen(str);i++)
    {
       str1[i-(j+1)]=str[i];        
    }
    str1[i]='\0';
    arith2=atoi(str1);
    arith1=arith1*pow((double)10,(double)arith2);
    return arith1;
}

void createNodeList(int node){
        if (root_node == NULL) {
            curr_node= (struct nodes*) malloc (sizeof(struct nodes));
            curr_node->name=node;
            curr_node->next=NULL;
            root_node=curr_node;
            return;
        }
        
        curr_node=root_node;
        while (curr_node != NULL){
            if (curr_node->name != node)
            {
                pre_node=curr_node;
                curr_node=curr_node->next;
            }
            else 
                break;             
        }    

        if (curr_node == NULL) {
            curr_node = (struct nodes*) malloc (sizeof(struct nodes));
            curr_node->name=node;
            curr_node->next=NULL;
            pre_node->next=curr_node;
        }
}

int count_nodes(){
    int counter = 0;
    curr_node=root_node;
    while (curr_node != NULL){
        counter++;
        curr_node=curr_node->next;
    }
    return counter;
}

void creatVoltList(FILE *fp){

	char *buf=(char*)malloc(sizeof(char)*buf_length);
        int node=0;
	VoltT *new;
	VoltT *curr;	

	new = (VoltT*) malloc(sizeof(VoltT));

	fscanf(fp,"%s",buf);
	new->name=(char*)malloc(sizeof(char)*(strlen(buf)+2));
        strcpy(new->name,"V");
	strcat(new->name, buf);
	
	fscanf(fp,"%d",&node);
	new->node1=node;
	if(node==0){
              ground=1;
        }
        createNodeList(node);
	
 	fscanf(fp,"%d",&node);
	new->node2=node;
	if(node==0){
            ground=1;
        }
        createNodeList(node);
	
	fscanf(fp,"%s",buf);
	if(strstr(buf, "e")==NULL) new->value=atof(buf);	
	else new->value=convertStringToDouble(buf);

        if(rootV==NULL){
		rootV=new;
		new->next=NULL;
	}
	else{
        	curr=rootV;
		while(curr->next!=NULL){
			curr=curr->next;
		}
		curr->next=new;
		new->next=NULL;
	}
	while((buf[0]=fgetc(fp))!='\n'&&(buf[0]!=EOF)){}
	sparse_elements++;
}

void creatAmberList(FILE *fp){

	char *buf=(char*)malloc(sizeof(char)*buf_length);
        int node=0;
	AmperT *new;
	AmperT *curr;

	new = (AmperT*) malloc(sizeof(AmperT));

	fscanf(fp,"%s",buf);
	new->name=(char*)malloc(sizeof(char)*(strlen(buf)+2));
        strcpy(new->name,"I");
	strcat(new->name, buf);
	
	fscanf(fp,"%d",&node);
	new->node1=node;
	if(node==0){
           ground=1;
        }
        createNodeList(node);
	
	fscanf(fp,"%d",&node);
	new->node2=node;
	if(node==0){
           ground=1;
        }
        createNodeList(node);
	
	fscanf(fp,"%s",buf);
	if(strstr(buf, "e")==NULL) new->value=atof(buf);	
	else new->value=convertStringToDouble(buf);

	if(rootI==NULL){
		rootI=new;
		new->next=NULL;
	}
	else{
		curr=rootI;
		while(curr->next!=NULL){
			curr=curr->next;
		}
		curr->next=new;
		new->next=NULL;
	}
	while((buf[0]=fgetc(fp))!='\n'&&(buf[0]!=EOF)){}
}

void creatResistanceList(FILE *fp){
	
	char *buf=(char*)malloc(sizeof(char)*buf_length);
        int node=0;
	ResistanceT *new;
	ResistanceT *curr;
	
	new = (ResistanceT*)malloc(sizeof(ResistanceT));
	
	fscanf(fp,"%s",buf);
	new->name=(char*)malloc(sizeof(char)*(strlen(buf)+2));
	strcpy(new->name,"R");
	strcat(new->name, buf);

	
	fscanf(fp,"%d",&node);
	new->node1=node;
	if(node==0){
           ground=1;
        }
        createNodeList(node);
	
	fscanf(fp,"%d",&node);
	new->node2=node;
	if(node==0){
           ground=1;
        }
        createNodeList(node);
	
	fscanf(fp,"%s",buf);
	if(strstr(buf, "e")==NULL) new->value=atof(buf);	
	else new->value=convertStringToDouble(buf);

	if(rootR==NULL){
		rootR=new;
		new->next=NULL;
	}
	else{
		curr=rootR;
		while(curr->next!=NULL){
			curr=curr->next;
		}

		curr->next=new;
		new->next=NULL;
        }
	while((buf[0]=fgetc(fp))!='\n'&&(buf[0]!=EOF)){}
	sparse_elements++;
}

void creatCapacitorList(FILE *fp){

	char *buf=(char*)malloc(sizeof(char)*buf_length);
        int node=0;
	CapacitorT *new;
	CapacitorT *curr;

	new = (CapacitorT*) malloc(sizeof(CapacitorT));

	fscanf(fp,"%s",buf);
	new->name=(char*)malloc(sizeof(char)*(strlen(buf)+2));
        strcpy(new->name,"C");
	strcat(new->name, buf);

	
	fscanf(fp,"%d",&node);
	new->node1=node;
	if(node==0){
           ground=1;
        }
        createNodeList(node);	

	fscanf(fp,"%d",&node);
	new->node2=node;
	if(node==0){
           ground=1;
        }
        createNodeList(node);
	
	fscanf(fp,"%s",buf);
	if(strstr(buf, "e")==NULL) new->value=atof(buf);	
	else new->value=convertStringToDouble(buf);

	if(rootC==NULL){
		rootC=new;
		new->next=NULL;
	}
	else{
		curr=rootC;
		while(curr->next!=NULL){
			curr=curr->next;
		}
		curr->next=new;
		new->next=NULL;
	}
	while((buf[0]=fgetc(fp))!='\n'&&(buf[0]!=EOF)){}
}

void creatInductorList(FILE *fp){

	char *buf=(char*)malloc(sizeof(char)*buf_length);
        int node=0;
	InductorT *new;
	InductorT *curr;	

	new = (InductorT*) malloc(sizeof(InductorT));

	fscanf(fp,"%s",buf);
	new->name=(char*)malloc(sizeof(char)*(strlen(buf)+2));
        strcpy(new->name,"L");
	strcat(new->name, buf);

	
	fscanf(fp,"%d",&node);
	new->node1=node;
	if(node==0){
           ground=1;
        }
        createNodeList(node);
	
	fscanf(fp,"%d",&node);
	new->node2=node;
	if(node==0){
           ground=1;
        }
        createNodeList(node);
	
	fscanf(fp,"%s",buf);
	if(strstr(buf, "e")==NULL) new->value=atof(buf);	
	else new->value=convertStringToDouble(buf);

	if(rootL==NULL){
		rootL=new;
		new->next=NULL;
	}
	else{
		curr=rootL;
		while(curr->next!=NULL){
			curr=curr->next;
		}

		curr->next=new;
		new->next=NULL;
	}
	while((buf[0]=fgetc(fp))!='\n'&&(buf[0]!=EOF)){}
	sparse_elements++;
}

void creatDiodeList(FILE *fp){
	
	char *buf=(char*)malloc(sizeof(char)*buf_length);
        int node=0;
	DiodeT *new;
	DiodeT *curr;

	new = (DiodeT*) malloc(sizeof(DiodeT));

	fscanf(fp,"%s",buf);
	new->name=(char*)malloc(sizeof(char)*(strlen(buf)+2));
        strcpy(new->name,"D");
	strcat(new->name, buf);

	
	fscanf(fp,"%d",&node);
	new->node1=node;
	if(node==0){
           ground=1;
        }
        createNodeList(node);
	
	fscanf(fp,"%d",&node);
	new->node2=node;
	if(node==0){
           ground=1;
        }
        createNodeList(node);
		
	if(rootD==NULL){
		rootD=new;
		new->next=NULL;
	}
	else{
		curr=rootD;
		while(curr->next!=NULL){
			curr=curr->next;
		}
		curr->next=new;
		new->next=NULL;
	}
	while((buf[0]=fgetc(fp))!='\n'&&(buf[0]!=EOF)){}

}

void creatMOSList(FILE *fp){

	char *buf=(char*)malloc(sizeof(char)*buf_length);
        int node=0;
	MosT *new;
	MosT *curr;

	new = (MosT*) malloc(sizeof(MosT));

	fscanf(fp,"%s",buf);
	new->name=(char*)malloc(sizeof(char)*(strlen(buf)+2));
        strcpy(new->name,"M");
	strcat(new->name, buf);

		
	fscanf(fp,"%d",&node);
	new->D=node;
	if(node==0){
           ground=1;
        }
        createNodeList(node);
	
	fscanf(fp,"%d",&node);
	new->G=node;
	if(node==0){
           ground=1;
        }
        createNodeList(node);
	
	fscanf(fp,"%d",&node);
	new->S=node;
	if(node==0){
           ground=1;
        }
        createNodeList(node);
	
	if(rootM==NULL){
		rootM=new;
		new->next=NULL;
	}
	else{
		curr=rootM;
		while(curr->next!=NULL){
			curr=curr->next;
		}
        	curr->next=new;
		new->next=NULL;
	}
        while((buf[0]=fgetc(fp))!='\n'&&(buf[0]!=EOF)){}
}

void creatBJTList(FILE *fp){

	char *buf=(char*)malloc(sizeof(char)*buf_length);
        int node=0;
	BjtT *new;
	BjtT *curr;

	new = (BjtT*) malloc(sizeof(BjtT));

	fscanf(fp,"%s",buf);
	new->name=(char*)malloc(sizeof(char)*(strlen(buf)+2));
        strcpy(new->name,"Q");
	strcat(new->name, buf);

	fscanf(fp,"%d",&node);
	new->C=node;
	if(node==0){
           ground=1;
        }
        createNodeList(node);
	
	fscanf(fp,"%d",&node);
	new->B=node;
	if(node==0){
           ground=1;
        }
        createNodeList(node);
	
	fscanf(fp,"%d",&node);
	new->E=node;
	if(node==0){
           ground=1;
        }
        createNodeList(node);
	
	if(rootB==NULL){
		rootB=new;
		new->next=NULL;
	}
	else{
		curr=rootB;
		while(curr->next!=NULL){
			curr=curr->next;
		}
		curr->next=new;
		new->next=NULL;
	}
	while((buf[0]=fgetc(fp))!='\n'&&(buf[0]!=EOF)){}
}

void printNodeList(){
    curr_node=root_node;
    while(curr_node!=NULL){
        printf("node: %d", curr_node->name);
        curr_node=curr_node->next;
    }
}

void printVoltList(VoltT *list){
   while(list!=NULL){
       printf(" %s  ",list->name);
       printf("node1 = %d  ",list->node1);
       printf("node2 = %d  ",list->node2);
       printf("value = %lf\n",list->value);
       list=list->next;
    }
}

void printAmperList(AmperT *list){
   while(list!=NULL){
       printf(" %s  ",list->name);
       printf("node1 = %d  ",list->node1);
       printf("node2 = %d  ",list->node2);
       printf("value = %lf\n",list->value);
       list=list->next;
    }
}

void printResistanceList(ResistanceT *list){
   while(list!=NULL){
       printf(" %s  ",list->name);
       printf("node1 = %d  ",list->node1);
       printf("node2 = %d  ",list->node2);
       printf("value = %lf\n",list->value);
       list=list->next;
    }
}

void printCapacitorList(CapacitorT *list){
   while(list!=NULL){
       printf(" %s  ",list->name);
       printf("node1 = %d  ",list->node1);
       printf("node2 = %d  ",list->node2);
       printf("value = %lf\n",list->value);
       list=list->next;
    }
}

void printInductorList(InductorT *list){
   while(list!=NULL){
       printf(" %s  ",list->name);
       printf("node1 = %d  ",list->node1);
       printf("node2 = %d  ",list->node2);
       printf("value = %lf\n",list->value);
       list=list->next;
    }
}

void printDiodeList(DiodeT *list){
   while(list!=NULL){
       printf(" %s  ",list->name);
       printf("node1 = %d  ",list->node1);
       printf("node2 = %d  ",list->node2);
       list=list->next;
    }
}

void printMosList(MosT *list){
   while(list!=NULL){
       printf(" %s  ",list->name);
       printf("drain = %d  ",list->D);
       printf("gate = %d ",list->G);
       printf("source = %d  ",list->S);
       list=list->next;
    }
}

void printBjttList(BjtT *list){
   while(list!=NULL){
       printf(" %s  ",list->name);
       printf("collector = %d  ",list->C);
       printf("base = %d  ",list->B);
       printf("emitter = %d  ",list->E);
       list=list->next;
    }
}

void analysis(FILE *fp){
    char *string;
    const char delimiters[] = " ,\t()\n";
    char line[1024];

    fgets(line,1024,fp);

    string= strtok (line, delimiters);

    if(strcmp(string,"OPTIONS")==0){
        //string = strtok (NULL, delimiters);
        while((string = strtok (NULL, delimiters))!=NULL){
            if ((string!=NULL) && (strstr(string,"SPD")!=NULL)){
               use_cholesky =1;
               use_lu=0;
            }
	    if (strstr(string,"SPARSE")!=NULL){
		sparse_option=1;
                printf("Using sparse\n");
                continue;
	    }
	    else if (strstr(string,"ITER")!=NULL){
               found_iter=1;
               itol_value=0.001;
               continue;
            }
	    else if (strstr(string,"ITOL")!=NULL){
               itol_value=atof(strtok (NULL, delimiters));
               continue;
            }
        }
        return;
    }
    else if(strcmp(string,"DC")==0){
        string = strtok (NULL, delimiters);
        if(string==NULL){
            return;
        }
        if(string[0]=='V'){
            int cnt;
            currV=rootV;
            cnt=count_nodes();
    
	        while(currV!=NULL){
                if(strcmp(currV->name,string+1)==0){
                    source=cnt;
                    break;
                }
                cnt++;
                currV=currV->next;  
            }
            start_value = atof(strtok (NULL, delimiters));
            end_value = atof(strtok (NULL, delimiters));
            step = atof(strtok (NULL, delimiters));
            found_dc_sweep=1;
            return;
        }
        else if(string[0]=='I')
		{//an meta to DC yparxei I, tote kanoume DC sweep allazontas kapoia pigi reumatos
      	    currI=rootI;

      	    while(currI!=NULL){
	    		if(!(strcmp(currI->name,(string+1)))){
	        		source=-1;	
	        		sweep_node1=currI->node1;	
	        		sweep_node2=currI->node2;	
	        		sweep_value=currI->value;
	        		break;
	    		}
	    		currI=currI->next;	
        	}
        	start_value = atof(strtok (NULL, delimiters));
     	 	end_value = atof(strtok (NULL, delimiters));
        	step = atof(strtok (NULL, delimiters));
        	found_dc_sweep=1;
        	return; 
        }
		else if(!(strcmp(string,"TRAN"))){
    	    return;
  		}
		else if(!(strcmp(string,"PLOT"))){
	    	return;
		}
    }
    else {
        return;
    }
}
