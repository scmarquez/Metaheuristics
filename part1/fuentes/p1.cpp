#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
/////////////////////////////////////////
//Generador de números pseudo-aleatorios a partir de una semilla Seed//
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>

#include <time.h>

using namespace std;

/***** GENERADOR DE NUMEROS PSEUDOALETORIOS *****/

/* used in random number generator below */
#define MASK 2147483647
#define PRIME 65539
#define SCALE 0.4656612875e-9

/*******************************************************************/
/*  Rand genera un numero real pseudoaleatorio entre 0 y 1,        */
/*  excluyendo el 1.						   */
/*  Randint genera un numero entero entre low y high, ambos 	   */
/*  incluidos.	   						   */
/*  Randfloat genera un numero real entre low y high, incluido low */
/*  y no incluido high                                             */
/*******************************************************************/


int Seed = 3;
#define Rand()  (( Seed = ( (Seed * PRIME) & MASK) ) * SCALE )

#define Randint(low,high) ( (int) (low + (high-(low)+1) * Rand()))

#define Randfloat(low,high) ( (low + (high-(low))*Rand()))

struct VectorEvaluation{
	vector<int>v;
	int evaluation;
	bool operator <(VectorEvaluation &ve);
	bool operator >(VectorEvaluation &ve);
};
bool VectorEvaluation:: operator <(VectorEvaluation &ve){
	if(evaluation < ve.evaluation){
		return true;
	}
	return false;
}
bool VectorEvaluation:: operator >(VectorEvaluation &ve){
	if(evaluation > ve.evaluation){
		return true;
	}
	return false;
}
struct PosEval{
	int pos;
	int eval;
	bool operator <(PosEval &pe);
	bool operator >(PosEval &pe);
};
bool PosEval:: operator <(PosEval &pe){
	if(eval < pe.eval){
		return true;
	}
	return false;
}
bool PosEval:: operator >(PosEval &pe){
	if(eval > pe.eval){
		return true;
	}
	return false;
}
struct FlowStruct{
	int unit;
	int flow;
	bool operator <(FlowStruct &fs);
	bool operator >(FlowStruct &fs);
};
struct DistanceStruct{
	int unit;
	int distance;
	bool operator <(DistanceStruct &ds);
	bool operator >(DistanceStruct &ds);
};
bool FlowStruct:: operator <(FlowStruct &fs){
	if(flow < fs.flow){
		return true;
	}
	else{
		return false;
	}
}
bool FlowStruct:: operator >(FlowStruct &fs){
	if(flow > fs.flow){
		return true;
	}
	else{
		return false;
	}
}
bool DistanceStruct:: operator <(DistanceStruct &ds){
	if(distance < ds.distance){
		return true;
	}
	else{
		return false;
	}
}
bool DistanceStruct:: operator >(DistanceStruct &ds){
	if(distance > ds.distance){
		return true;
	}
	else{
		return false;
	}
}
class Data{
	private: 
		int size;
		vector<vector<int>> distances;
		vector<vector<int>> flows;
	public:
		Data(){}
		void LoadFromFile(string filename){
			ifstream entrada;
			entrada.open(filename);
			if(entrada.good()){
				entrada >> size;
				distances.resize(size);
				flows.resize(size);
				for(int i = 0; i < size; i++){
					distances[i].resize(size);
					flows[i].resize(size);
				}
				for(int i = 0;i < size; i++){
					for(int j = 0; j < size; j++){
						entrada >>distances[i][j];
					}
				}
				for(int i  =0; i < size; i++){
					for(int j =0; j < size; j++){
						entrada >>flows[i][j];
					}
				}
				entrada.close();
			}
			else{
				cerr << "Can not open such file, please close all the programs that are using it and retry.\n";
			}
			
		}
		vector<vector<int>> GetDistances(){
			return distances;
		}
		vector<vector<int>> GetFlows(){
			return flows;
		}
		int GetFlowIn(int i, int j){
			if(i < size && j < size){
				return flows[i][j];
			}
			else{
				cerr << "No aviable position\n";
			}
		}
		int GetDistanceIn(int i, int j){
			if(i < size && j < size){
				return distances[i][j];
			}
			else{
				cerr << "No aviable position\n";
			}
		}
		int GetSize(){
			return size;
		}
};
void GenerateRandomVector(vector<int> &v){
	
	int size = v.size();
	for(int i = 0; i < size; i++){
		v[i] = i;
	}
	//Barajar
	for(int i = 0; i < size; i++){
		int p1 = Randint(i,size-1);

		int swap = v[i];
		v[i] = v[p1];
		v[p1] = swap;
	}
}
vector<int> greedy(Data &data){
	vector<FlowStruct>flows; 		//Vector que mantiene la influencia de cada elemento
	vector<DistanceStruct>distances; //Vector que mantiene la centralidad de cada elemento
	vector<int>solution;
	//Calculo de la influencia de cada elemento(matriz de flujo)
	int size = data.GetSize();
	flows.resize(size);
	distances.resize(size);
	solution.resize(size);

	for(int i = 0; i < data.GetSize(); i++){
		//Para el elemento i-ésimo
		int val = 0;
		for(int j = 0; j < data.GetSize(); j++){
			val += data.GetFlowIn(i,j);
		}
		FlowStruct fs;
		fs.unit = i;
		fs.flow = val;
		flows[i] = fs;
	}
	sort(flows.begin(),flows.end());
	reverse(flows.begin(),flows.end());
	//Cálculo de las distancias entre secciones
	for(int i = 0; i < data.GetSize(); i++){
		//Para la posición i-ésima
		int val = 0;
		for(int j = 0; j < data.GetSize(); j++){
			val+=data.GetDistanceIn(i,j);
		}
		DistanceStruct ds;
		ds.unit = i;
		ds.distance = val;
		distances[i] = ds;
	}
	sort(distances.begin(),distances.end());
	//construlle la solución
	for(int i = 0; i < size; i++){
		solution[distances[i].unit] = flows[i].unit;
	}
	return solution;
}
int EvaluateVector(vector<int> &v, Data &data){
	//v contiene en cada posición la localización asignada a cada unidad
	int eval = 0;
	int size = data.GetSize();
	for(int i = 0; i < size; i++){
		for(int j = 0; j < size; j++){
			eval += data.GetFlowIn(i,j)*data.GetDistanceIn(v[i],v[j]);
		}
	}
	return eval;
}
int EvaluateNeigbour(int eval, int &p1, int &p2,Data &data,vector<int> &v){
	//Hay que invertir el cambio
	int swap = v[p1];
	v[p1] = v[p2];
	v[p2] = swap;


	//k = i, s = p1, r = p2
	int incremento = 0;
	int size = data.GetSize();
	for(int i = 0; i < size; i++){
		if(i != p1 && i != p2){
			incremento += data.GetFlowIn(p2,i)*(data.GetDistanceIn(v[p1],v[i])-data.GetDistanceIn(v[p2],v[i]))+
			data.GetFlowIn(p1,i)*(data.GetDistanceIn(v[p2],v[i])-data.GetDistanceIn(v[p1],v[i]))+
			data.GetFlowIn(i,p2)*(data.GetDistanceIn(v[i],v[p1])-data.GetDistanceIn(v[i],v[p2]))+
			data.GetFlowIn(i,p1)*(data.GetDistanceIn(v[i],v[p2])-data.GetDistanceIn(v[i],v[p1]));
		}
	}
	//Deshacer los cambios en el vector

	swap = v[p1];
	v[p1] = v[p2];
	v[p2] = swap;

	eval += incremento;
	return eval;
}
bool FindNextBestNeighbour(vector<int> &v, Data &data,vector<bool> &DLB,int &iterations, int maxIterations, int &ev){
	//Evalua inicialmente el vector de partida y obtiene la mejor solución actual(CurrentSol)
	int CurrenrSol = EvaluateVector(v,data);
	//Obtiene eltamaño del problema
	int size = data.GetSize();
	//Indica que no ha encontrado ninguna solución aún
	bool found = false;

	int swap;
	//Usa un vector aleatorio para saber que orden usa para las permutaciones
	vector<int>order;
	order.resize(size);
	GenerateRandomVector(order);

	for(int i = 0; i < size && !found; i++){
		//Obtiene la primera posición que usará para las permutaciones()posi
		int posi = order[i];
		//Comprueba la máscara de bits para saber si hace la permutación con otro elemento
		if(DLB[posi]){
			for(int j = 0; j < size && !found; j++){
				if(iterations < maxIterations){
					//Obtiene el elemento con el que se permutará posi(posj)
					int posj = order[j];
					if(i != j){
						//Intercambia posi con posj en el vector y lo evalua
						swap = v[posi];
						v[posi] = v[posj];
						v[posj] = swap;
						int eval = EvaluateNeigbour(CurrenrSol,posi,posj,data,v);
						//int eval = EvaluateVector(v,data);
						//Incrementa el número de llamadas a la evaluación
						iterations++;

						if(eval < CurrenrSol){
							//Si la evaluación es positiva, pone el found a true y sale del bucle
							CurrenrSol = eval;
							found = true;
							//Resitua el valor de máscara a positivos
							DLB[posj] = true;
							DLB[posi] = true;
						}
						else{
							//Deja el vector como estaba si la evaluacion es negativa.
							swap = v[posi];
							v[posi] = v[posj];
							v[posj] = swap;
						}
					}
				}			
			}	
		}
		if(!found){
			//Si después de salir del bucle de contador j tiene el found a false, se pone la máscara a false.
			DLB[posi] = false;
		}
		
	}
	ev = CurrenrSol;
	return found;
}
vector<int> LocalSearch(Data &data, int nIter){
	//Parámetros básicos del problema
	int size = data.GetSize();
	bool Improve = true;
	int CurrentIters = 0;
	int CurrenrSol;
	//Es necesario crear el vector de partida de forma aleatoria
	vector<int> InitSolution;
	InitSolution.resize(size);
	GenerateRandomVector(InitSolution);
	//Máscara Dont lool bits inicializada a true(puede mirar todas las posiciones)
	vector<bool> DLB;
	DLB.resize(size);
	for(int i = 0; i < size; i++){
		DLB[i] = true;
	}
	//Mientras encuentre una mejor solución y no supere el número de iteraciones busca más vecinos
	while(CurrentIters < nIter && Improve){
		/////////////////////////////////////////////////////////////
		int ev;
		Improve = FindNextBestNeighbour(InitSolution,data,DLB,CurrentIters,nIter,ev);
		////////////////////////////////////////////////////////////
	}
	return InitSolution;
}
void PositionCross(vector<int> &parent1, vector<int> &parent2,vector<int> &son){
	int size = parent1.size();
	//Busca las posiciones comunes de los padres y las copia en hijo.
	//Para saber que valores no han sido usados se declara una máscara de uso
	vector<bool>used;
	used.resize(size);
	//Inicializa el vector a -1 para controlar las posiciones sin usar
	for(int i = 0; i < size; i++){
		son[i] = -1;
	}
	for(int i = 0; i < size; i++){
		int p1 = parent1[i];
		int p2 = parent2[i]; 
		//Si ambas posiciones son iguales se copia en el hijo(son[i])
		if(p1 == p2){
			son[i] = p1;
			used[p1] = true;
		}
		else{
			used[p1] = false;
		}
	}
	//Una vez que se tienen todas las posiciones comunes se completa con los restantes de forma aleatoria.
	vector<int>pos;
	pos.resize(size);
	GenerateRandomVector(pos);
	int write_pos = 0;
	for(int i = 0; i < size; i++){
		if(!used[pos[i]]){
			//Busca la siguinte posicion de escritura
			while(son[write_pos] != -1){
				write_pos++;
			}
			son[write_pos] = pos[i];
		}
	}
} 
void OXCross(vector<int> &parent1, vector<int> &parent2 ,vector<int> &son, float porcentaje){
	//Escoge la subcadena central hasta que alcanza el porcentaje especificado
	int size = parent1.size();
	int n_elements = size*porcentaje;
	int halfPos = size/2;
	int inicio,fin;
	
	inicio = halfPos - n_elements/2;
	fin = halfPos + n_elements/2;

	vector<bool> escogidos;
	escogidos.resize(size);
	//Inicializa la máscara de elementos. 
	for(int i = 0; i < size; i++){
		if(i >= inicio && i <= fin){
			escogidos[parent1[i]] = true;
		}
		else{
			escogidos[parent1[i]] = false;
		}
	}
	//Una vez inicializada la máscara se rellena el vector hijo
	int last_not_used = 0;
	for(int i = 0; i < size; i++){
		if(i >= inicio && i <= fin){
			son[i] = parent1[i];
		}
		else{
			//Busca el siguiente elemento de parent2 con la máscara en false
			while(escogidos[parent2[last_not_used]]){
				last_not_used++;
			}
			son[i] = parent2[last_not_used];
			last_not_used++;
		}
		
	}
}
vector<int> GeneticoGeneracional(Data &data,int population_size,int n_iters, float prob_mutation, float prob_son,int n_parents,bool ox){
	int size = data.GetSize();
	int iters = 0;
	//Se crea la población inicial
	//La población inicial está formada por 2 vectores, uno de soluciones(vector de vectorEvalution)
	//y otro de indices con evaluaciones(vector de PosEval)
	vector<VectorEvaluation> population;
	population.resize(population_size);

	vector<PosEval> population_index;
	population_index.resize(population_size);
	//Se rellenan los vectores con soluciones aleatorias
	for(int i = 0; i < population_size; i++){
		//Se inicializa el vector de soluciones
		population[i].v.resize(size);
		GenerateRandomVector(population[i].v);
		//Se añade su índice y su evaluación
		int ev = EvaluateVector(population[i].v,data);
		population_index[i].pos = i;
		population_index[i].eval = ev;
		population[i].evaluation = ev;
	}
	
	while(iters < n_iters){
		///////////////////////////////////////////////////////////////
		//////////SELECCIÓN DE PADRES//////////////////////////////////
		///////////////////////////////////////////////////////////////
		//Para seleccionar a los padres se crean 2 pares de vectores de poseval con índices aleatoios(candidatos1,candidatos2)
		//Se hace un torneo binario entre cada posición i del vector y el ganador va al vector de padres.
		vector<int>candidatos1;
		candidatos1.resize(n_parents);
		for(int i = 0; i < n_parents; i++){
			candidatos1[i] = Randint(0,population_size-1);
		}
		vector<int>candidatos2;
		candidatos2.resize(n_parents);
		for(int i = 0; i < n_parents; i++){
			candidatos2[i] = Randint(0,population_size-1);
		}
		//Ahora se aplica un torneo binario para pareja de padres y se inserta en la lista de padres
		vector<int>padres;
		padres.resize(n_parents);
		for(int i = 0; i < n_parents; i++){
			if(population[candidatos1[i]].evaluation < population[candidatos2[i]].evaluation){
				//Si gana el elemento de candidatos 1 se introduce ese índice en el vector
				padres[i] = candidatos1[i];
			}
			else{
				padres[i] = candidatos2[i];
			}
		}
		//cout << "Genera padres\n";
		//Una vez calculados los padres se crean los vectores asociados a los hijos
		vector<VectorEvaluation>sons_population;
		vector<PosEval>sons_index;
		//Primero se calcula el número estimado de hijos que se van a obtener
		int n_sons = ceil((n_parents/2)*prob_son);
		sons_population.resize(n_sons);
		sons_index.resize(n_sons);
		for(int i = 0; i < n_sons;i++){
			sons_population[i].v.resize(size);
		}
		//Ahora se cruzan los padres para obtener un hijo
		for(int i = 0; i < n_sons && iters < n_iters; i++){
			if(!ox){
				PositionCross(population[padres[i*2]].v, population[padres[i*2+1]].v, sons_population[i].v);
			}
			else{
				OXCross(population[padres[i*2]].v, population[padres[i*2+1]].v, sons_population[i].v,0.45);
			}
			
			sons_population[i].evaluation = EvaluateVector(sons_population[i].v,data);
			iters++;
			if(iters == n_iters){
				n_sons = i;
			}
			//Ahora se rellena el vector de índice correpondiente
			sons_index[i].pos = i;
			sons_index[i].eval = sons_population[i].evaluation;
		}
		//cout << "Genera hijos\n";
		//Mutación de los descendientes.
		//Se calcula el número estimado de mutaciones
		
		int n_mutations = n_sons*size*prob_mutation;
		//cout << "Espera hacer " << n_mutations << " mutaciones.\n";
		for(int i = 0; i < n_mutations; i++){
			
			int individuo = Randint(0,n_sons-1);
			//cout << "individuo: " << individuo << endl;
			
			int gen1 = Randint(0,size-1);
			//cout << "gen1: " << gen1 << endl;
			
			int gen2 =Randint(0,size-1);
			//cout << "gen2: " << gen2 << endl;
			//Solución antes de mutar al individuo.
			int CurrenrSol = sons_population[individuo].evaluation;

			int swap = sons_population[individuo].v[gen1];
			sons_population[individuo].v[gen1] = sons_population[individuo].v[gen2];
			sons_population[individuo].v[gen2] = swap;

			//Se reevalua el elemento mutado
			int ev = EvaluateNeigbour(CurrenrSol,gen1,gen2,data,sons_population[individuo].v);
			//int ev = EvaluateVector(sons_population[individuo].v,data);
			sons_population[individuo].evaluation = ev;
			sons_index[individuo].eval =  ev;
			
		}
		//cout << "Genera mutaciones\n";
		
		//Ahora se tiene que crear la nueva población. 
		//Es necesario ordenar los vectores
		sort(population_index.begin(), population_index.end());
		sort(sons_index.begin(), sons_index.end());

		int pos_insertado = population_size - n_sons;
		//Se escogen los mejores hijos hasta que no queden hijos o hasta que se rellene por completo la población.
		if(pos_insertado > 0){
			//En este caso se tendrá que completar con la población anterior, por lo que el mejor anterior siempre sobrevive.
			for(int i = 0; i < n_sons && i < population_size; i++){
				//Los hijos se inseertarán el el vector de población 
				int indice_population = population_index[pos_insertado+i].pos;
				int indice_sons = sons_index[i].pos;
				population[indice_population] = sons_population[indice_sons];
			}
		}
		else{
			//En este caso la población de hijos es mayor o igual que la anterior.
			//Se copian los mejores hijos y en último lugar si el mejor de la generación anterior
			//es distinto del mejor de los hijos también se añade.
			//Esta comprobación previa puede ahorrar la siguiente más costosa
			if(population_index[0].eval < sons_index[0].pos || population_index[0].eval > sons_index[0].pos){
				//En este caso ya sabemos que es distinto y simplemente copiamos todos los hijos respetando
				//al mejor padre
				for(int i = 1; i < population_size; i++){
					population[population_index[i].pos] = sons_population[sons_index[i-1].pos];
				}
			}
			else{
				//Ahora si es necesario comprobar si la solución es igual
				bool equal = true;
				for(int i = 0; i < size && equal; i++){
					int element1 = population[population_index[0].pos].v[i];
					int element2 = sons_population[sons_index[0].pos].v[i];
					if(element2 != element1){
						equal = false;
					}
				}
				if(equal){
					//no se respeta al mejor anterior, porque ya está
					for(int i = 0; i < population_size; i++){
						population[i] = sons_population[sons_index[i-1].pos];
					}
				}
				else{
					//se respeta al mejor anterior
					for(int i = 1; i < population_size; i++){
						population[population_index[i].pos] = sons_population[sons_index[i-1].pos];
					}
				}
			}
		}
		//Reinicializar el vector de indices
		for(int i = 0; i < population_size; i++){
			population_index[i].pos = i;
			population_index[i].eval = population[i].evaluation;
		}
		//cout << "crea la nueva población y termina la iteración " << iters << endl;
	}

	sort(population_index.begin(),population_index.end());
	return population[population_index[0].pos].v;
}
vector<int>GeneticoEstacionario(Data &data, int population_size,int n_iters,float prob_mutation,bool ox){
	int size = data.GetSize();
	int iters = 0;
	//Se crea la población inicial
	//La población inicial está formada por 2 vectores, uno de soluciones(vector de vectorEvalution)
	vector<VectorEvaluation> population;
	population.resize(population_size);

	//Se rellenan los vectores con soluciones aleatorias
	for(int i = 0; i < population_size; i++){
		//Se inicializa el vector de soluciones
		population[i].v.resize(size);
		GenerateRandomVector(population[i].v);
		//Se añade su índice y su evaluación
		int ev = EvaluateVector(population[i].v,data);
		population[i].evaluation = ev;
	}
	//Para realizar las mutaciones se calcula la frecuencia de mutación estimada(inverso de la probabilidad)
	int frecuencia_mutacion = 1/prob_mutation;
	int mutaciones = 0;
	while(iters < n_iters){
		//Se seleccionan las dos parejas de padres que serán cruzadas
		//Se aplican 2 torneos binarios por pareja

		//Se escogen 4 candidatos al azar para la primera pareja de padres
		int candidato1 = Randint(0,population_size-1);
		int candidato2 = Randint(0,population_size-1);
		int candidato3 = Randint(0,population_size-1);
		int candidato4 = Randint(0,population_size-1);

		int padre1;
		int padre2;
		//Torneo binario 
		if(population[candidato1].evaluation < population[candidato2].evaluation){
			padre1 = candidato1;
		}
		else{
			padre1 = candidato2;
		}
		if(population[candidato3].evaluation < population[candidato4].evaluation){
			padre2 = candidato3;
		}
		else{
			padre2 = candidato4;
		}

		//Ahora se crea la segunda pareja
		candidato1 = Randint(0,population_size-1);
		candidato2 = Randint(0,population_size-1);
		candidato3 = Randint(0,population_size-1);
		candidato4 = Randint(0,population_size-1);

		int padre3;
		int padre4;
		//Torneo binario 
		if(population[candidato1].evaluation < population[candidato2].evaluation){
			padre3 = candidato1;
		}
		else{
			padre3 = candidato2;
		}
		if(population[candidato3].evaluation < population[candidato4].evaluation){
			padre4 = candidato3;
		}
		else{
			padre4 = candidato4;
		}

		//Una vez escogidos los padres se realizan los cruces.
		//Primero se crean los vectores necesarios para los hijos
		VectorEvaluation hijo1;
		VectorEvaluation hijo2;
		hijo1.v.resize(size);
		hijo2.v.resize(size);
		if(ox){
			OXCross(population[padre1].v,population[padre2].v,hijo1.v,0.45);
			OXCross(population[padre3].v,population[padre4].v,hijo2.v,0.45);
		}
		else{
			PositionCross(population[padre1].v,population[padre2].v,hijo1.v);
			PositionCross(population[padre3].v,population[padre4].v,hijo2.v);
		}
		//Ahora se mutan los hijos 
		//Generar un hijo supone crear "size" genes nuevos susceptibles a la mutación
		mutaciones += 2* size;
		if(mutaciones > frecuencia_mutacion){
			//Ahora deberíamos mutar algún elemento de los hijos
			int hijo = Randint(0,1);
			int gen1 = Randint(0,size-1);
			int gen2 = Randint(0,size-1);

			if(hijo == 0){
				//Se muta el hijo1
				int swap = hijo1.v[gen1];
				hijo1.v[gen1] = hijo1.v[gen2];
				hijo1.v[gen2] = swap;
			}
			else{
				//Se muta el hijo2
				int swap = hijo2.v[gen1];
				hijo2.v[gen1] = hijo2.v[gen2];
				hijo2.v[gen2] = swap;
			}
		}
		//Se evaluan los dos hijos generados
		hijo1.evaluation = EvaluateVector(hijo1.v,data);
		hijo2.evaluation = EvaluateVector(hijo2.v,data);
		iters +=2;
		//Una vez se tiene los hijos creados se evaluan y se hace un torneo con los peores de la población
		//Búsqueda de los dos peores de la población
		int indice_peor = 0;
		int indice_segundo_peor = 0;

		for(int i = 0; i < population_size; i++){
			if(population[i].evaluation > population[indice_peor].evaluation){
				indice_peor = i;
			}
			if(population[i].evaluation > population[indice_segundo_peor].evaluation && indice_peor != i){
				indice_segundo_peor = i;
			}
		}

		vector<PosEval> torneo;
		torneo.resize(4);
		torneo[0].pos = indice_segundo_peor;
		torneo[0].eval = population[indice_segundo_peor].evaluation;
		torneo[1].pos = indice_peor;
		torneo[1].eval = population[indice_peor].evaluation;
		torneo[2].eval = hijo1.evaluation;
		torneo[2].pos = -1;
		torneo[3].eval = hijo2.evaluation;
		torneo[3].pos = -2;

		sort(torneo.begin(),torneo.end());

		if(torneo[0].pos == -1){
			//mete al hijo1 donde el peor
			population[indice_peor] = hijo1;
			if(torneo[1].pos == -2){
				//mete al hijo 2 donde el segundo peor
				population[indice_segundo_peor] = hijo2;
			}
		}
		else{
			if(torneo[0].pos == -2){
				//mete al hijo2 donde el peor
				population[indice_peor] = hijo2;
				if(torneo[1].pos == -1){
					//mete al hijo 1 donde el segundo peor
					population[indice_segundo_peor] = hijo1;
				}
			}
		}
	}
	//Antes de terminar busca el mejor resultado y lo devuelve
	int indice_mejor = 0;
	for(int i = 0; i < population_size; i++){
		if(population[i].evaluation < population[indice_mejor].evaluation){
			indice_mejor = i;
		}
	}
	return population[indice_mejor].v;
}
void LocalSearchInicializada(Data &data, int nIter, vector<int> &v, int &cIters,int &eval){
	//Parámetros básicos del problema
	int size = data.GetSize();
	bool Improve = true;
	int CurrentIters = 0;
	int CurrenrSol;
	
	//Máscara Dont lool bits inicializada a true(puede mirar todas las posiciones)
	vector<bool> DLB;
	DLB.resize(size);
	for(int i = 0; i < size; i++){
		DLB[i] = true;
	}
	//Mientras encuentre una mejor solución y no supere el número de iteraciones busca más vecinos
	while(CurrentIters < nIter && Improve){
		/////////////////////////////////////////////////////////////
		Improve = FindNextBestNeighbour(v,data,DLB,CurrentIters,nIter,eval);
		////////////////////////////////////////////////////////////
	}
	
	cIters += CurrentIters;
	//cout << CurrentIters << endl;
}
vector<int>Memetic1(Data &data,int population_size, int n_iters, float prob_mutation, float prob_son,int n_parents,bool ox){
	int size = data.GetSize();
	int iters = 0;
	int generations = 0;
	//Se crea la población inicial
	//La población inicial está formada por 2 vectores, uno de soluciones(vector de vectorEvalution)
	//y otro de indices con evaluaciones(vector de PosEval)
	vector<VectorEvaluation> population;
	population.resize(population_size);

	vector<PosEval> population_index;
	population_index.resize(population_size);
	//Se rellenan los vectores con soluciones aleatorias
	for(int i = 0; i < population_size; i++){
		//Se inicializa el vector de soluciones
		population[i].v.resize(size);
		GenerateRandomVector(population[i].v);
		//Se añade su índice y su evaluación
		int ev = EvaluateVector(population[i].v,data);
		population_index[i].pos = i;
		population_index[i].eval = ev;
		population[i].evaluation = ev;
	}
	
	while(iters < n_iters){
		///////////////////////////////////////////////////////////////
		//////////SELECCIÓN DE PADRES//////////////////////////////////
		///////////////////////////////////////////////////////////////
		//Para seleccionar a los padres se crean 2 pares de vectores de poseval con índices aleatoios(candidatos1,candidatos2)
		//Se hace un torneo binario entre cada posición i del vector y el ganador va al vector de padres.
		vector<int>candidatos1;
		candidatos1.resize(n_parents);
		for(int i = 0; i < n_parents; i++){
			candidatos1[i] = Randint(0,population_size-1);
		}
		vector<int>candidatos2;
		candidatos2.resize(n_parents);
		for(int i = 0; i < n_parents; i++){
			candidatos2[i] = Randint(0,population_size-1);
		}
		//Ahora se aplica un torneo binario para pareja de padres y se inserta en la lista de padres
		vector<int>padres;
		padres.resize(n_parents);
		for(int i = 0; i < n_parents; i++){
			if(population[candidatos1[i]].evaluation < population[candidatos2[i]].evaluation){
				//Si gana el elemento de candidatos 1 se introduce ese índice en el vector
				padres[i] = candidatos1[i];
			}
			else{
				padres[i] = candidatos2[i];
			}
		}
		//cout << "Genera padres\n";
		//Una vez calculados los padres se crean los vectores asociados a los hijos
		vector<VectorEvaluation>sons_population;
		vector<PosEval>sons_index;
		//Primero se calcula el número estimado de hijos que se van a obtener
		int n_sons = ceil((n_parents/2)*prob_son);
		sons_population.resize(n_sons);
		sons_index.resize(n_sons);
		for(int i = 0; i < n_sons;i++){
			sons_population[i].v.resize(size);
		}
		//Ahora se cruzan los padres para obtener un hijo
		for(int i = 0; i < n_sons && iters < n_iters; i++){
			if(!ox){
				PositionCross(population[padres[i*2]].v, population[padres[i*2+1]].v, sons_population[i].v);
			}
			else{
				OXCross(population[padres[i*2]].v, population[padres[i*2+1]].v, sons_population[i].v,0.45);
			}
			
			sons_population[i].evaluation = EvaluateVector(sons_population[i].v,data);
			iters++;
			if(iters == n_iters){
				n_sons = i;
			}
			//Ahora se rellena el vector de índice correpondiente
			sons_index[i].pos = i;
			sons_index[i].eval = sons_population[i].evaluation;
		}
		//cout << "Genera hijos\n";
		//Mutación de los descendientes.
		//Se calcula el número estimado de mutaciones
		
		int n_mutations = n_sons*size*prob_mutation;
		//cout << "Espera hacer " << n_mutations << " mutaciones.\n";
		for(int i = 0; i < n_mutations; i++){
			
			int individuo = Randint(0,n_sons-1);
			//cout << "individuo: " << individuo << endl;
			
			int gen1 = Randint(0,size-1);
			//cout << "gen1: " << gen1 << endl;
			
			int gen2 =Randint(0,size-1);
			//cout << "gen2: " << gen2 << endl;
			//Solución antes de mutar al individuo.
			int CurrenrSol = sons_population[individuo].evaluation;

			int swap = sons_population[individuo].v[gen1];
			sons_population[individuo].v[gen1] = sons_population[individuo].v[gen2];
			sons_population[individuo].v[gen2] = swap;

			//Se reevalua el elemento mutado
			int ev = EvaluateNeigbour(CurrenrSol,gen1,gen2,data,sons_population[individuo].v);
			//int ev = EvaluateVector(sons_population[individuo].v,data);
			sons_population[individuo].evaluation = ev;
			sons_index[individuo].eval =  ev;
			
		}
		//cout << "Genera mutaciones\n";
		
		//Ahora se tiene que crear la nueva población. 
		//Es necesario ordenar los vectores
		sort(population_index.begin(), population_index.end());
		sort(sons_index.begin(), sons_index.end());

		int pos_insertado = population_size - n_sons;
		//Se escogen los mejores hijos hasta que no queden hijos o hasta que se rellene por completo la población.
		if(pos_insertado > 0){
			//En este caso se tendrá que completar con la población anterior, por lo que el mejor anterior siempre sobrevive.
			for(int i = 0; i < n_sons && i < population_size; i++){
				//Los hijos se inseertarán el el vector de población 
				int indice_population = population_index[pos_insertado+i].pos;
				int indice_sons = sons_index[i].pos;
				population[indice_population] = sons_population[indice_sons];
			}
		}
		
		else{
			//En este caso la población de hijos es mayor o igual que la anterior.
			//Se copian los mejores hijos y en último lugar si el mejor de la generación anterior
			//es distinto del mejor de los hijos también se añade.
			//Esta comprobación previa puede ahorrar la siguiente más costosa
			if(population_index[0].eval < sons_index[0].pos || population_index[0].eval > sons_index[0].pos){
				//En este caso ya sabemos que es distinto y simplemente copiamos todos los hijos respetando
				//al mejor padre
				for(int i = 1; i < population_size; i++){
					population[population_index[i].pos] = sons_population[sons_index[i-1].pos];
				}
			}
			else{
				//Ahora si es necesario comprobar si la solución es igual
				bool equal = true;
				for(int i = 0; i < size && equal; i++){
					int element1 = population[population_index[0].pos].v[i];
					int element2 = sons_population[sons_index[0].pos].v[i];
					if(element2 != element1){
						equal = false;
					}
				}
				if(equal){
					//no se respeta al mejor anterior, porque ya está
					for(int i = 0; i < population_size; i++){
						population[i] = sons_population[sons_index[i-1].pos];
					}
				}
				else{
					//se respeta al mejor anterior
					for(int i = 1; i < population_size; i++){
						population[population_index[i].pos] = sons_population[sons_index[i-1].pos];
					}
				}
			}
		}

		generations++;
		if(generations == 10){
			//Reemplazar toda la población por la bl
			for(int i = 0; i < population_size; i++){
				int ev;
				LocalSearchInicializada(data,400,population[i].v,iters,ev);
				population[i].evaluation = ev;
			}
			//Reinicializar el vector de indices
			for(int i = 0; i < population_size; i++){
				population_index[i].pos = i;
				population_index[i].eval = population[i].evaluation;
			}
			generations = generations%10;
		}
	}
	sort(population.begin(),population.end());
	return population[0].v;
}
vector<int>Memetic2(Data &data,int population_size, int n_iters, float prob_mutation, float prob_son,int n_parents,float prob_bl,bool ox){
	int size = data.GetSize();
	int iters = 0;
	int generations = 0;
	//Se crea la población inicial
	//La población inicial está formada por 2 vectores, uno de soluciones(vector de vectorEvalution)
	//y otro de indices con evaluaciones(vector de PosEval)
	vector<VectorEvaluation> population;
	population.resize(population_size);

	vector<PosEval> population_index;
	population_index.resize(population_size);
	//Se rellenan los vectores con soluciones aleatorias
	for(int i = 0; i < population_size; i++){
		//Se inicializa el vector de soluciones
		population[i].v.resize(size);
		GenerateRandomVector(population[i].v);
		//Se añade su índice y su evaluación
		int ev = EvaluateVector(population[i].v,data);
		population_index[i].pos = i;
		population_index[i].eval = ev;
		population[i].evaluation = ev;
	}
	
	while(iters < n_iters){
		///////////////////////////////////////////////////////////////
		//////////SELECCIÓN DE PADRES//////////////////////////////////
		///////////////////////////////////////////////////////////////
		//Para seleccionar a los padres se crean 2 pares de vectores de poseval con índices aleatoios(candidatos1,candidatos2)
		//Se hace un torneo binario entre cada posición i del vector y el ganador va al vector de padres.
		vector<int>candidatos1;
		candidatos1.resize(n_parents);
		for(int i = 0; i < n_parents; i++){
			candidatos1[i] = Randint(0,population_size-1);
		}
		vector<int>candidatos2;
		candidatos2.resize(n_parents);
		for(int i = 0; i < n_parents; i++){
			candidatos2[i] = Randint(0,population_size-1);
		}
		//Ahora se aplica un torneo binario para pareja de padres y se inserta en la lista de padres
		vector<int>padres;
		padres.resize(n_parents);
		for(int i = 0; i < n_parents; i++){
			if(population[candidatos1[i]].evaluation < population[candidatos2[i]].evaluation){
				//Si gana el elemento de candidatos 1 se introduce ese índice en el vector
				padres[i] = candidatos1[i];
			}
			else{
				padres[i] = candidatos2[i];
			}
		}
		//cout << "Genera padres\n";
		//Una vez calculados los padres se crean los vectores asociados a los hijos
		vector<VectorEvaluation>sons_population;
		vector<PosEval>sons_index;
		//Primero se calcula el número estimado de hijos que se van a obtener
		int n_sons = ceil((n_parents/2)*prob_son);
		sons_population.resize(n_sons);
		sons_index.resize(n_sons);
		for(int i = 0; i < n_sons;i++){
			sons_population[i].v.resize(size);
		}
		//Ahora se cruzan los padres para obtener un hijo
		for(int i = 0; i < n_sons && iters < n_iters; i++){
			if(!ox){
				PositionCross(population[padres[i*2]].v, population[padres[i*2+1]].v, sons_population[i].v);
			}
			else{
				OXCross(population[padres[i*2]].v, population[padres[i*2+1]].v, sons_population[i].v,0.45);
			}
			
			sons_population[i].evaluation = EvaluateVector(sons_population[i].v,data);
			iters++;
			if(iters == n_iters){
				n_sons = i;
			}
			//Ahora se rellena el vector de índice correpondiente
			sons_index[i].pos = i;
			sons_index[i].eval = sons_population[i].evaluation;
		}
		//cout << "Genera hijos\n";
		//Mutación de los descendientes.
		//Se calcula el número estimado de mutaciones
		
		int n_mutations = n_sons*size*prob_mutation;
		//cout << "Espera hacer " << n_mutations << " mutaciones.\n";
		for(int i = 0; i < n_mutations; i++){
			
			int individuo = Randint(0,n_sons-1);
			//cout << "individuo: " << individuo << endl;
			
			int gen1 = Randint(0,size-1);
			//cout << "gen1: " << gen1 << endl;
			
			int gen2 =Randint(0,size-1);
			//cout << "gen2: " << gen2 << endl;
			//Solución antes de mutar al individuo.
			int CurrenrSol = sons_population[individuo].evaluation;

			int swap = sons_population[individuo].v[gen1];
			sons_population[individuo].v[gen1] = sons_population[individuo].v[gen2];
			sons_population[individuo].v[gen2] = swap;

			//Se reevalua el elemento mutado
			int ev = EvaluateNeigbour(CurrenrSol,gen1,gen2,data,sons_population[individuo].v);
			//int ev = EvaluateVector(sons_population[individuo].v,data);
			sons_population[individuo].evaluation = ev;
			sons_index[individuo].eval =  ev;
			
		}
		//cout << "Genera mutaciones\n";
		
		//Ahora se tiene que crear la nueva población. 
		//Es necesario ordenar los vectores
		sort(population_index.begin(), population_index.end());
		sort(sons_index.begin(), sons_index.end());

		int pos_insertado = population_size - n_sons;
		//Se escogen los mejores hijos hasta que no queden hijos o hasta que se rellene por completo la población.
		if(pos_insertado > 0){
			//En este caso se tendrá que completar con la población anterior, por lo que el mejor anterior siempre sobrevive.
			for(int i = 0; i < n_sons && i < population_size; i++){
				//Los hijos se inseertarán el el vector de población 
				int indice_population = population_index[pos_insertado+i].pos;
				int indice_sons = sons_index[i].pos;
				population[indice_population] = sons_population[indice_sons];
			}
		}
		
		else{
			//En este caso la población de hijos es mayor o igual que la anterior.
			//Se copian los mejores hijos y en último lugar si el mejor de la generación anterior
			//es distinto del mejor de los hijos también se añade.
			//Esta comprobación previa puede ahorrar la siguiente más costosa
			if(population_index[0].eval < sons_index[0].pos || population_index[0].eval > sons_index[0].pos){
				//En este caso ya sabemos que es distinto y simplemente copiamos todos los hijos respetando
				//al mejor padre
				for(int i = 1; i < population_size; i++){
					population[population_index[i].pos] = sons_population[sons_index[i-1].pos];
				}
			}
			else{
				//Ahora si es necesario comprobar si la solución es igual
				bool equal = true;
				for(int i = 0; i < size && equal; i++){
					int element1 = population[population_index[0].pos].v[i];
					int element2 = sons_population[sons_index[0].pos].v[i];
					if(element2 != element1){
						equal = false;
					}
				}
				if(equal){
					//no se respeta al mejor anterior, porque ya está
					for(int i = 0; i < population_size; i++){
						population[i] = sons_population[sons_index[i-1].pos];
					}
				}
				else{
					//se respeta al mejor anterior
					for(int i = 1; i < population_size; i++){
						population[population_index[i].pos] = sons_population[sons_index[i-1].pos];
					}
				}
			}
		}
		generations++;
		if(generations == 10){
			//Seleccionar sobre que individuos se lanza la bl
			int n_busquedas_locales = population_size*prob_bl;

			for(int i = 0; i < n_busquedas_locales; i++){
				int indice_bl = Randint(0,population_size-1);
				int ev;
				LocalSearchInicializada(data,400,population[indice_bl].v,iters,ev);
				population[indice_bl].evaluation = ev;
			}
			//Reinicializar el vector de indices
			for(int i = 0; i < population_size; i++){
				population_index[i].pos = i;
				population_index[i].eval = population[i].evaluation;
			}
			generations = generations % 10;
		}
		
	}
	sort(population.begin(),population.end());
	return population[0].v;
}

vector<int>Memetic3(Data &data,int population_size, int n_iters, float prob_mutation, float prob_son,int n_parents,float prob_bl,bool ox){
	int size = data.GetSize();
	int iters = 0;
	int generations = 0;
	//Se crea la población inicial
	//La población inicial está formada por 2 vectores, uno de soluciones(vector de vectorEvalution)
	//y otro de indices con evaluaciones(vector de PosEval)
	vector<VectorEvaluation> population;
	population.resize(population_size);

	vector<PosEval> population_index;
	population_index.resize(population_size);
	//Se rellenan los vectores con soluciones aleatorias
	for(int i = 0; i < population_size; i++){
		//Se inicializa el vector de soluciones
		population[i].v.resize(size);
		GenerateRandomVector(population[i].v);
		//Se añade su índice y su evaluación
		int ev = EvaluateVector(population[i].v,data);
		population_index[i].pos = i;
		population_index[i].eval = ev;
		population[i].evaluation = ev;
	}
	
	while(iters < n_iters){
		///////////////////////////////////////////////////////////////
		//////////SELECCIÓN DE PADRES//////////////////////////////////
		///////////////////////////////////////////////////////////////
		//Para seleccionar a los padres se crean 2 pares de vectores de poseval con índices aleatoios(candidatos1,candidatos2)
		//Se hace un torneo binario entre cada posición i del vector y el ganador va al vector de padres.
		vector<int>candidatos1;
		candidatos1.resize(n_parents);
		for(int i = 0; i < n_parents; i++){
			candidatos1[i] = Randint(0,population_size-1);
		}
		vector<int>candidatos2;
		candidatos2.resize(n_parents);
		for(int i = 0; i < n_parents; i++){
			candidatos2[i] = Randint(0,population_size-1);
		}
		//Ahora se aplica un torneo binario para pareja de padres y se inserta en la lista de padres
		vector<int>padres;
		padres.resize(n_parents);
		for(int i = 0; i < n_parents; i++){
			if(population[candidatos1[i]].evaluation < population[candidatos2[i]].evaluation){
				//Si gana el elemento de candidatos 1 se introduce ese índice en el vector
				padres[i] = candidatos1[i];
			}
			else{
				padres[i] = candidatos2[i];
			}
		}
		//cout << "Genera padres\n";
		//Una vez calculados los padres se crean los vectores asociados a los hijos
		vector<VectorEvaluation>sons_population;
		vector<PosEval>sons_index;
		//Primero se calcula el número estimado de hijos que se van a obtener
		int n_sons = ceil((n_parents/2)*prob_son);
		sons_population.resize(n_sons);
		sons_index.resize(n_sons);
		for(int i = 0; i < n_sons;i++){
			sons_population[i].v.resize(size);
		}
		//Ahora se cruzan los padres para obtener un hijo
		for(int i = 0; i < n_sons && iters < n_iters; i++){
			if(!ox){
				PositionCross(population[padres[i*2]].v, population[padres[i*2+1]].v, sons_population[i].v);
			}
			else{
				OXCross(population[padres[i*2]].v, population[padres[i*2+1]].v, sons_population[i].v,0.45);
			}
			
			sons_population[i].evaluation = EvaluateVector(sons_population[i].v,data);
			iters++;
			if(iters == n_iters){
				n_sons = i;
			}
			//Ahora se rellena el vector de índice correpondiente
			sons_index[i].pos = i;
			sons_index[i].eval = sons_population[i].evaluation;
		}
		//cout << "Genera hijos\n";
		//Mutación de los descendientes.
		//Se calcula el número estimado de mutaciones
		
		int n_mutations = n_sons*size*prob_mutation;
		//cout << "Espera hacer " << n_mutations << " mutaciones.\n";
		for(int i = 0; i < n_mutations; i++){
			
			int individuo = Randint(0,n_sons-1);
			//cout << "individuo: " << individuo << endl;
			
			int gen1 = Randint(0,size-1);
			//cout << "gen1: " << gen1 << endl;
			
			int gen2 =Randint(0,size-1);
			//cout << "gen2: " << gen2 << endl;
			//Solución antes de mutar al individuo.
			int CurrenrSol = sons_population[individuo].evaluation;

			int swap = sons_population[individuo].v[gen1];
			sons_population[individuo].v[gen1] = sons_population[individuo].v[gen2];
			sons_population[individuo].v[gen2] = swap;

			//Se reevalua el elemento mutado
			int ev = EvaluateNeigbour(CurrenrSol,gen1,gen2,data,sons_population[individuo].v);
			//int ev = EvaluateVector(sons_population[individuo].v,data);
			sons_population[individuo].evaluation = ev;
			sons_index[individuo].eval =  ev;
			
		}
		//cout << "Genera mutaciones\n";
		
		//Ahora se tiene que crear la nueva población. 
		//Es necesario ordenar los vectores
		sort(population_index.begin(), population_index.end());
		sort(sons_index.begin(), sons_index.end());

		int pos_insertado = population_size - n_sons;
		//Se escogen los mejores hijos hasta que no queden hijos o hasta que se rellene por completo la población.
		if(pos_insertado > 0){
			//En este caso se tendrá que completar con la población anterior, por lo que el mejor anterior siempre sobrevive.
			for(int i = 0; i < n_sons && i < population_size; i++){
				//Los hijos se inseertarán el el vector de población 
				int indice_population = population_index[pos_insertado+i].pos;
				int indice_sons = sons_index[i].pos;
				population[indice_population] = sons_population[indice_sons];
			}
		}
		
		else{
			//En este caso la población de hijos es mayor o igual que la anterior.
			//Se copian los mejores hijos y en último lugar si el mejor de la generación anterior
			//es distinto del mejor de los hijos también se añade.
			//Esta comprobación previa puede ahorrar la siguiente más costosa
			if(population_index[0].eval < sons_index[0].pos || population_index[0].eval > sons_index[0].pos){
				//En este caso ya sabemos que es distinto y simplemente copiamos todos los hijos respetando
				//al mejor padre
				for(int i = 1; i < population_size; i++){
					population[population_index[i].pos] = sons_population[sons_index[i-1].pos];
				}
			}
			else{
				//Ahora si es necesario comprobar si la solución es igual
				bool equal = true;
				for(int i = 0; i < size && equal; i++){
					int element1 = population[population_index[0].pos].v[i];
					int element2 = sons_population[sons_index[0].pos].v[i];
					if(element2 != element1){
						equal = false;
					}
				}
				if(equal){
					//no se respeta al mejor anterior, porque ya está
					for(int i = 0; i < population_size; i++){
						population[i] = sons_population[sons_index[i-1].pos];
					}
				}
				else{
					//se respeta al mejor anterior
					for(int i = 1; i < population_size; i++){
						population[population_index[i].pos] = sons_population[sons_index[i-1].pos];
					}
				}
			}
		}
		generations++;
		if(generations == 10){
			//Seleccionar sobre que individuos se lanza la bl
			int n_busquedas_locales = population_size*prob_bl;

			//Reinicializar el vector de indices
			for(int i = 0; i < population_size; i++){
				population_index[i].pos = i;
				population_index[i].eval = population[i].evaluation;
			}
			sort(population_index.begin(),population_index.end());
			for(int i = 0; i < n_busquedas_locales; i++){
				int ev;
				LocalSearchInicializada(data,400,population[population_index[i].pos].v,iters,ev);
				population[population_index[i].pos].evaluation = ev;
			}
			//Reinicializar el vector de indices
			for(int i = 0; i < population_size; i++){
				population_index[i].pos = i;
				population_index[i].eval = population[i].evaluation;
			}
			generations = generations % 10;
		}
		
	}
	sort(population.begin(),population.end());
	return population[0].v;
}


int main(int argc, char *argv[]){
	if(argc != 3){
		cerr << "Se esperaba como argumento el nombre del archivo y la semilla\n";
	}
	else{
		Seed = atoi(argv[2]);

		Data d;

		d.LoadFromFile(argv[1]);
		vector<int> v1;
		
		cout << argv[1] << " ";
		clock_t inicio;
		clock_t fin;
		
		inicio = clock();

		v1 = LocalSearch(d,50000);
		fin = clock();
		cout << EvaluateVector(v1,d)<< " ";
		cout << ((float)fin - float(inicio))/CLOCKS_PER_SEC << endl;
		
		inicio = clock();
		v1 = GeneticoGeneracional(d,50,50000,0.001,0.7,50,false);
		fin = clock();
		cout << EvaluateVector(v1,d)<< " ";
		cout << ((float)fin - float(inicio))/CLOCKS_PER_SEC << endl;
		

		
		inicio = clock();
		v1 = GeneticoGeneracional(d,50,50000,0.001,0.7,50,true);
		fin = clock();
		cout << EvaluateVector(v1,d)<<" ";
		cout << ((float)fin - float(inicio))/CLOCKS_PER_SEC << "\n";
		
		
		inicio = clock();
		v1 = GeneticoEstacionario(d,50,50000,0.001,false);
		fin = clock();
		cout << EvaluateVector(v1,d)<<" ";
		cout << ((float)fin - float(inicio))/CLOCKS_PER_SEC << endl;
		
		
		inicio = clock();
		v1 = GeneticoEstacionario(d,50,50000,0.001,true);
		fin = clock();
		cout << EvaluateVector(v1,d)<<" ";
		cout << ((float)fin - float(inicio))/CLOCKS_PER_SEC << endl;
		
		inicio = clock();
		v1 = Memetic1(d,50,50000,0.001,0.7,50,false);
		fin = clock();
		cout << EvaluateVector(v1,d)<<" ";
		cout << ((float)fin - float(inicio))/CLOCKS_PER_SEC << endl;
		
		
		inicio = clock();
		v1 = Memetic2(d,50,50000,0.001,0.7,50,0.1,false);
		fin = clock();
		cout << EvaluateVector(v1,d)<<" ";
		cout << ((float)fin - float(inicio))/CLOCKS_PER_SEC << endl;
		
		
		inicio = clock();
		v1 = Memetic3(d,50,50000,0.001,0.7,50,0.1,false);
		fin = clock();
		cout << EvaluateVector(v1,d)<<" ";
		cout << ((float)fin - float(inicio))/CLOCKS_PER_SEC << endl;
		
		inicio = clock();

		v1 = greedy(d);
		fin = clock();
		cout << EvaluateVector(v1,d)<< " ";
		cout << ((float)fin - float(inicio))/CLOCKS_PER_SEC << endl;
	}
}