#include<bits/stdc++.h>
#include<chrono>

using namespace std;
using namespace std::chrono;

std::random_device r;
std::seed_seq seed{r()};
std::mt19937 gen(seed);


//ofstream file("maxclus4.dat");
// variables for network building
int m; //number of hands for each node
int m0; //initial number of seeds, all connected to each other.
int N0;//final number of nodes
int N; // N =N0-m0
int bond_number;
int filenum;
int mm; //number of bonds to compare
int ensemble_count;

vector<int> ptr;
vector<int>bond_id;
vector<int>sbond_id;

vector < vector <int> > Nodes_vect; //collection of nodes, it will contain N elements, N= number of nodes
vector <int> NodeA; //node at one side of a bond, it will contain M elements, M= number of bonds
vector <int> NodeB; //node at another side of a bond, it will contain M elements, M= number of bonds

string largestclus1;


void seed_initialization()
{
	
	Nodes_vect.resize(N0);
	
	// Create seed nodes and connect them with each other
    for (int i=0; i<m0; i++) 
	{
        for (int j= i+1; j<m0; j++) 
		{
           Nodes_vect[i].push_back(j);
            Nodes_vect[j].push_back(i);
        } 
    }
	
	//cout<<"seed ok"<<endl;
	
	//seed bonds initialization
	for(int i=0; i<m0-1; i++)
	{
		for(int j=i+1; j<m0; j++)
		{
			NodeA.push_back(i);
			NodeB.push_back(j);
		}
	}
}



void network_build(double p)
{
	
	// Randomly select a node from the existing nodes and connect it to m neighbours
    for (int i = m0; i <N0; i++) 
	{
        // Select a random node from the existing nodes
        uniform_int_distribution<int> distribution(0, i-1);
        int selected_node = distribution(gen);
        
        uniform_real_distribution<double> random(0,1);
        double random_number = random(gen);
        
        // Select m neighbours of the selected node
        vector<int> neighbours = Nodes_vect[selected_node];
        //shuffle(neighbours.begin(), neighbours.end(),gen);
       // neighbours.resize(m);
        Nodes_vect[i].push_back(selected_node);
        Nodes_vect[selected_node].push_back(i);
        NodeA.push_back(i);
		NodeB.push_back(selected_node);
        // Connect the new node to the selected m neighbours
        for (int j = 0; j < neighbours.size(); j++) {
        	uniform_real_distribution<long double> random(0,1);
        double random_number = random(gen);
        
        	if (random_number <= p ){
			
            Nodes_vect[i].push_back(neighbours[j]);
            Nodes_vect[neighbours[j]].push_back(i);
			NodeA.push_back(i);
			NodeB.push_back(neighbours[j]);
        }
		}
    }
}





int findroot(int i)
{
    if (ptr[i]<0) return i;
    return ptr[i] = findroot(ptr[i]);
}



void ensemble_initialization()
{
	 for(int a=0; a<NodeA.size(); a++)
  {
      bond_id.push_back(a);
  }
	
    for(int h=0; h<bond_id.size(); h++)
    {
        sbond_id.push_back(bond_id[h]);
    }
	shuffle(sbond_id.begin(),sbond_id.end(), gen);

//cout<<"Ensemble initialization ok"<<endl;
}

vector<long double>product_of_clustersize;


int bondselection(int A)
{


   for(int i=0; i<mm; i++)
   {
      int e1 = sbond_id[A+ i];
      int s1 = NodeA[e1];
      int s2 = NodeB[e1];
      int s1_root = findroot(s1);
      int s2_root = findroot(s2);
	  int s1_size= -ptr[s1_root];
      int s2_size = -ptr[s2_root];
	  long double product = (long double) s1_size* (long double) s2_size;
	 
	  
	  
      product_of_clustersize.push_back(product);
   }

  
  long double selected_bond_size = product_of_clustersize[0];
   int selected_bond = A;
   for(int j=1; j<mm; j++)
   {
       if(product_of_clustersize[j]<selected_bond_size)
       {
          selected_bond_size= product_of_clustersize[j];
           selected_bond= A+j;

       }
   }


    int original_selection = sbond_id[selected_bond];

    swap(sbond_id[A],sbond_id[selected_bond]);
    uniform_int_distribution<int> dist(A+1,sbond_id.size()-1);

    for(int k=1; k<mm; k++)
        {
            int f = dist(gen);
            swap(sbond_id[A+k], sbond_id[f]);
        }

    vector<long double>().swap(product_of_clustersize);
	
	

    return original_selection;

}

void percolation()
{   
     vector<long double> big(N0, 0);
    ofstream file(largestclus1);
    for (int ensemble = 0; ensemble < ensemble_count; ensemble++)
    {  // cout<<"ok1"<<endl;
        seed_initialization();
        network_build(0.8);
        //cout<<"ok2"<<endl;
        ensemble_initialization();
        //cout<<"ok3"<<endl;
        for (int i=0; i<N0; i++)
	    {
		ptr.push_back(-1);
	    }
	    //ptr.resize(N0,-1);
        int bond_number = sbond_id.size();
         
        int current_big = 0;
        for (int a = 0; a < N0; a++)
        {
            int x, y, a1;

            if (a <= bond_number - mm)
            {
                a1 = bondselection(a);
            }
            else
            {
                a1 = sbond_id[a];
            }

            x = NodeA[a1];
            y = NodeB[a1];

            int x1 = findroot(x);
            int y1 = findroot(y);

            if (ptr[x1] == -1 && ptr[y1] == -1)
            {
                ptr[x1] += ptr[y1];
                ptr[y1] = x1;
            }
            else if (x1 != y1)
            {
                if (ptr[x1] > ptr[y1])
                {
                    ptr[y1] += ptr[x1];
                    ptr[x1] = y1;
                }
                else
                {
                    ptr[x1] += ptr[y1];
                    ptr[y1] = x1;
                }
            }

            if (-ptr[findroot(x)] > current_big)
            {
                current_big = -ptr[findroot(x)];
            }
             big[a] +=(double)current_big/N0;
             
              
             //cout<<"ok4"<<endl;
        }

         

        // Reset the network and bond configurations for each ensemble
       /* NodeA.clear();
        NodeB.clear();
        ptr.clear();
        bond_id.clear();
        sbond_id.clear();*/
        
        vector<int>().swap(NodeA);
        vector<int>().swap(NodeB);
        vector<int>().swap(ptr);
        vector<int>().swap(bond_id);
        vector<int>().swap(sbond_id);
        vector<vector<int> >().swap(Nodes_vect);
        
		cout<<"ensemble no."<<ensemble<<" "<<"is done"<<endl;    
    }
   
    
    // Calculate average cluster sizes and print the results
    for (int i = 0; i < N0; i++)
    {
         double avg_cluster_size = (double)(big[i]) / (ensemble_count)  ;
        
        //cout<<fixed<< setprecision(10) << (double)(i) / N0 << "  " << setprecision(10) << avg_cluster_size << endl;
        file <<fixed<<setprecision(10) << (double)(i) / N0 << "  " <<setprecision(10) << avg_cluster_size << endl;
       
    }
  
   
}



int main(int argc, char* argv[]){
	mm=atoi(argv[1]);
	N0 =atoi(argv[2]);
    m0=atoi(argv[3]);
    ensemble_count=atoi(argv[4]);
    filenum=atoi(argv[5]);
    
    largestclus1="mda_mm"+to_string(mm)+"_N0"+to_string(N0)+"_m0"+to_string(m0)+"_en"+to_string(ensemble_count)+"_filenum"+to_string(filenum)+".dat";
    
	percolation();
	
    return 0;
}

