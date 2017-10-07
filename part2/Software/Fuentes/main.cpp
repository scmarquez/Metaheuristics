#include <iostream>
#include <vector>
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

using namespace std;
struct UnitCostStruct{
    int unit;
    int cost;
    bool operator <(UnitCostStruct &ucs);
    bool operator >(UnitCostStruct &ucs);
};
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

bool UnitCostStruct:: operator <(UnitCostStruct &ucs){
    if(cost < ucs.cost){
        return true;
    }
    else{
        return false;
    }
}
bool UnitCostStruct:: operator >(UnitCostStruct &ucs){
    if(cost > ucs.cost){
        return true;
    }
    else{
        return false;
    }
}
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
        //solution[distances[i].unit] = flows[i].unit;
        solution[flows[i].unit] = distances[i].unit;
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
/************************************************/
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
    //Máscara Dont look bits inicializada a true(puede mirar todas las posiciones)
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
/************************************************/
void BarajarVectorParejas(vector<pair<int,int>> &v){
    int size = v.size();
    for(int i= 0; i < size; i++){
        int p1 = Randint(i,size-1);
        pair<int,int> swap;
        swap = v[i];
        v[i] = v[p1];
        v[p1] = swap;
    }
}

void GenerateRandomPair(int size, vector<pair<int,int>> &v){
    v.resize((size*size-size)/2);
    int pos  =0;
    for(int i = 0;i < size; i++){
        for (int j = i+1; j < size; j++) {
            v[pos].first = i;
            v[pos].second = j;
            pos++;
        }
    }
    BarajarVectorParejas(v);
}

vector<int> SimulatedAnnealing(Data &d, int max_evals){
    int size = d.GetSize();

    int max_vecinos = 10*size;//Número máximo de vecinos que se pueden generar
    int max_exitos = floor(0.1*max_vecinos);//Número máximo de soluciones que se pueden aceptar
    vector<int> mejor_solucion; //Mejor solución obtenida hasta el momento
    vector<int> solucion_actual;//Solución de partida del problema y con la que se trabaja
    solucion_actual.resize(size);
    GenerateRandomVector(solucion_actual);//Primera solución de la que partimos
    //Mejor solución será igual a la primera genrada
    mejor_solucion = solucion_actual;
    int mejor_eval = EvaluateVector(mejor_solucion,d);
    bool exito = true; //Si en un enfriamiento(bucle interno) no se encuentra nada mejor, se pone a false
    int exitos = 0;
    int vecinos_generados = 0;
    int n_evals = 0;
    int evaluacion_actual = mejor_eval;

    //Esquema de enfriamiento
    float mu = 0.4;
    float fi = 0.4;
    float t0 = mu*evaluacion_actual/-log(fi);
    float tf = pow(10,-3);
    int M = (float)max_evals/(float)max_vecinos;//Enfriamientos a realizar
    float B;//Variable necesaria para el enfriamiento de Cauchy
    B = (t0-tf)/(M*t0*tf);

    float t = t0; //Temperatura del sistema
    int iteracion = 1; //Contador de iteraciones actuales
    int n_entradas_por_aleatorio = 0;


    while(n_evals < max_evals && exito){
        //Mientras el algoritmo mejore y el número de evaluaciones sea menor al establecido
        exitos = 0;//Cuenta el número de vecinos generados en cada enfriamiento
        for(int i = 0; i < max_vecinos && exitos < max_exitos;i++) {
            //Se genera de forma aleatoria el intercambio
            int p1;
            int p2;
            p1 = Randint(0,size-1);
            p2 = Randint(0,size-1);
            //para la pareja i-ésima
            int Swap = solucion_actual[p1];
            solucion_actual[p1] = solucion_actual[p2];
            solucion_actual[p2] = Swap;
            //Se comprueba la evaluación de la nueva solución creada
            int evaluacion_nueva = EvaluateNeigbour(evaluacion_actual, p1, p2, d, solucion_actual);
            n_evals++;
            //Ahora se decide si se acepta o si se rechaza la nueva solución
            if (evaluacion_nueva < evaluacion_actual) {
                //Si la nueva solución es mejor que la actual, se acepta sin más
                exitos++;
                //Ahora comprobamos si es la mejor solución encontrada
                if(evaluacion_nueva < mejor_eval){
                    mejor_solucion = solucion_actual;
                    mejor_eval = evaluacion_nueva;
                }
                evaluacion_actual = evaluacion_nueva;
            }
            else {
                //En este caso se aceptará una peor solución si la temperatura es lo bastante alta
                if (exp((-1 * (evaluacion_nueva - evaluacion_actual)) / (iteracion * t)) >= Randfloat(0, 1)) {

                    n_entradas_por_aleatorio++;

                    //Entonces puede escoger la nueva solución

                    exitos++;
                    if(evaluacion_nueva < mejor_eval){
                        mejor_solucion = solucion_actual;
                        mejor_eval = evaluacion_actual;
                    }
                    evaluacion_actual = evaluacion_nueva;
                }
                else{
                    //No escoge la solución actual
                    //Es necesario invertir los cambios realizados en el vector
                    Swap = solucion_actual[p1];
                    solucion_actual[p1] = solucion_actual[p2];
                    solucion_actual[p2] = Swap;
                }
            }
            vecinos_generados++;
        }
        if(exitos == 0){
            exito = false;
        }
        iteracion++;
        //Genera una nueva temperatura para el sistema

        t = t/(1+(B*t));
        //t = 0.99*t;

    }
    return mejor_solucion;
}

vector<int>RandomGreedy(Data &data,float umbral_calidad){

    vector<FlowStruct>flows; 		 //Vector que mantiene la influencia de cada elemento
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
    /*
    //Para cada elemento se crea una lista con los mejores candidatos posibles y se escoge uno de forma aleatoria
    //Es necesario contemplar las localizaciones asignadas
    vector<bool> usedLocations;
    usedLocations.resize(size);
    for(int i = 0; i < size; i++){
        usedLocations[i] = false;
    }
    vector<bool> usedUnits;
    usedUnits.resize(size);
    for(int i = 0; i < size; i++){
        usedUnits[i] = false;
    }


    //Se añaden 2 asignaciones al vector de soluciones
    //Se añade la mejor unidad a la mejor localización y la segunda mejor a la segunda mejor localización
    solution[distances[0].unit] = flows[0].unit;
    usedLocations[distances[0].unit] = true;
    solution[distances[1].unit] = flows[1].unit;
    usedLocations[distances[1].unit] = true;
    int current_eval = data.GetFlowIn(flows[0].unit,flows[1].unit)*data.GetDistanceIn(distances[0].unit,distances[1].unit);
    //Ahora se asigan el resto de unidades a cada una de las localizaciones y se elige una aleatoria dentro de un umbral
    for(int i = 2; i < size; i++){
        int n_candidatos = 0;
        vector<int>vector_evaluaciones;
        //El tamaño de este vector será igual al número de localizaciones libres(size-i)
        vector_evaluaciones.resize(size);
        for(int j = 0; j < size; j++){
            vector_evaluaciones[j] = 0;
        }

        //Para la localización i-ésima en orden
        //Se prueba a asignarle cada unidad posible
        int eval = current_eval;
        for(int k = 2; k < size; k++){
            //Si la localización no ha sido usada puede probarse
            if(usedUnits[distances[k].unit] == false){
                //Se recalcula la evaluación si se asignase esta unidad
                for(int l = 0;l < size; l++){
                    if(usedLocations[l] == true){
                        eval += data.GetDistanceIn(l,i)*data.GetFlowIn(solution[l],k);
                    }
                }
                //Se añade el calculo al vector de evaluaciones
                vector_evaluaciones[k] = eval;
                n_candidatos++;
            }
            //Se calcula el umbral(busca el mejor y el peor de los costes)
            int indice_mejor = 0;
            int indice_peor = 0;
            for(int j = 0; j < size; j++){
                if(vector_evaluaciones[j] != -1 && vector_evaluaciones[j] > vector_evaluaciones[indice_peor]){
                    //Este es el peor encontrado hasta el momento
                    indice_peor = j;
                }
                else if(vector_evaluaciones[j] != -1 && vector_evaluaciones[j] < vector_evaluaciones[indice_peor]){
                    //Este es el mejor encontrado hasta el momento
                    indice_mejor = j;
                }
            }
            int umbral = vector_evaluaciones[indice_mejor] + umbral_calidad*(vector_evaluaciones[indice_peor]-vector_evaluaciones[indice_mejor]);
            //Se crea el vector de candidatos
            int indice_candidatos = 0;
            vector<int> candidatos;
            candidatos.resize(n_candidatos);
            for(int j = 0; j < size; j++){
                if(vector_evaluaciones[j] != -1 && vector_evaluaciones[j] < umbral) {
                    //En este caso se acepta como candidato y se añade j al vector de candidatos
                    candidatos[indice_candidatos] = j;
                    indice_candidatos++;
                }
            }
            //Se escoge una de forma aleatoria
            int unidad = candidatos[Randint(0,n_candidatos-1)];
            solution[distances[i].unit] = unidad;
            usedLocations[i] = true;
        }
    }
    */
    //Se inicializa el vecor solución
    for(int i = 0;i < size; i++){
        solution[i] = -1;
    }
    //Se escogen los 2 primeros elementos
    solution[flows[0].unit] = distances[0].unit;
    solution[flows[1].unit] = distances[1].unit;
    //Es necesario llevar un registro de las localizaciones ya asignadas
    vector<bool> used_locations;//registro de localizaciones ya asignadas a una localización
    used_locations.resize(size);
    for(int i = 0; i < size; i++){
        used_locations[i] = false;
    }
    used_locations[distances[0].unit] = true;
    used_locations[distances[1].unit] = true;

    int n_used_locations = 2;
    //Una vez que estas dos unidades están asignadas se calcula el coste de la solución parcial
    int cost = data.GetFlowIn(flows[0].unit,flows[1].unit)*data.GetDistanceIn(distances[0].unit,distances[1].unit);
    //Ahora por cada unidad restante se le asigna una unidad aleatoria dentro de un umbral de calidad
    for(int i = 2; i < size; i++){
        //Para la unidad i-esima(en orden)
        int unidad = flows[i].unit;
        //Se buscan las posibles localizaciones que puedan asignarse a esta unidad sin estar por encima del umbral
        vector<UnitCostStruct>localizaciones_candidatas; //par formado por first->unidad second->coste de asignacion de esa unidad
        localizaciones_candidatas.resize(size-n_used_locations);
        int localizaciones_index = 0;//Índice usado para insertar elementos en unidades_candidatas
        for(int j = 0; j < size; j++){
            int localizacion = distances[j].unit;
            //Si la unidad no se ha usado la podemos tener en cuenta
            if(!used_locations[localizacion]){
                //Se calcula el coste de la unidad
                int unit_cost = cost;
                for(int k = 0; k < size; k++){
                    //Si esa unidad tiene una localización asiganda
                    int s = solution[k];
                    if(s != -1){
                        int distance = data.GetDistanceIn(localizacion,solution[k]);
                        int flow = data.GetFlowIn(unidad,k);
                        unit_cost += distance*flow;
                    }
                }
                //Se asigna la localización al vector de candidatos

                UnitCostStruct ucs;
                ucs.cost = unit_cost;
                ucs.unit = localizacion;
                localizaciones_candidatas[localizaciones_index] = ucs;
                localizaciones_index++;
            }
        }
        sort(localizaciones_candidatas.begin(),localizaciones_candidatas.end());
        //Calculo del umbral
        int coste_mejor = localizaciones_candidatas[0].cost;
        int coste_peor = localizaciones_candidatas[size-n_used_locations-1].cost;
        float umbral = coste_mejor + umbral_calidad*(coste_peor-coste_mejor);
        //Se establece el índice de corte
        int indice_corte = 0;
        bool found = false;
        for(int j = 0; j < size-n_used_locations && !found; j++){
            if(localizaciones_candidatas[j].cost > umbral){
                indice_corte = j-1;
                found = true;
            }
        }
        //Una vez se encuentra el índice de corte se asigna una unidad entre la posición 0 y el índice de corte de forma aleatoria
        int asignacion = Randint(0,indice_corte);
        solution[unidad] = localizaciones_candidatas[asignacion].unit;
        used_locations[localizaciones_candidatas[asignacion].unit] = true;
        n_used_locations++;
        cost += localizaciones_candidatas[asignacion].cost;
    }
    return solution;
}

void LocalSearchInicializada(Data &data, int nIter, vector<int> &v, int &cIters,int &eval){
    //Parámetros básicos del problema
    int size = data.GetSize();
    bool Improve = true;
    int CurrentIters = 0;

    //Máscara Dont look bits inicializada a true(puede mirar todas las posiciones)
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

vector<int>Grasp(Data &d, int max_evals, int n_grasp){
    vector<int> solucion_actual;
    solucion_actual.resize(d.GetSize());
    solucion_actual = RandomGreedy(d,0.3);
    int evaluacion_actual = EvaluateVector(solucion_actual,d);
    vector<int> mejor_solucion = solucion_actual;
    int mejor_eval = evaluacion_actual;
    for(int i = 0; i < n_grasp; i++){
        int c_iters = 0;
        LocalSearchInicializada(d,50000,solucion_actual,c_iters,evaluacion_actual);
        evaluacion_actual = EvaluateVector(solucion_actual,d);
        if(evaluacion_actual < mejor_eval){
            mejor_solucion = solucion_actual;
            mejor_eval = evaluacion_actual;
        }
        solucion_actual = RandomGreedy(d,0.3);
    }
    return mejor_solucion;
}

vector<int>ILS(Data &d, int max_evals, int n_iters){
    int size = d.GetSize();
    int size_sublista = size/4;

    vector<int> solucion_actual;
    solucion_actual.resize(size);
    GenerateRandomVector(solucion_actual);
    int evaluacion_actual = EvaluateVector(solucion_actual,d);
    int current_evals = 0;
    LocalSearchInicializada(d,max_evals,solucion_actual,current_evals,evaluacion_actual);
    vector<int> mejor_solucion = solucion_actual;
    int mejor_eval = evaluacion_actual;


    for(int i = 1; i < n_iters;i++){
        //Se muta la solucion actual y se vuelve a lanzar la busqueda local.
        //Proceso de mutación
        //Primero se escoge el inicio de la sublista
        int inicio_sublista = Randint(0,size-size_sublista-1);

        //Se crea el vector de posiciones que serán alteradas a partir de la posición de inicio
        vector<int>posiciones;
        posiciones.resize(size_sublista);
        GenerateRandomVector(posiciones);
        for(int j = 0; j < size_sublista; j++){
            int swap = solucion_actual[j+inicio_sublista];
            solucion_actual[j+inicio_sublista] = solucion_actual[posiciones[j]+inicio_sublista];
            solucion_actual[posiciones[j]+inicio_sublista] = swap;
        }
        //BL sobre la solución mutada
        evaluacion_actual = EvaluateVector(solucion_actual,d);
        current_evals = 0;
        LocalSearchInicializada(d,max_evals,solucion_actual,current_evals,evaluacion_actual);
        //evaluacion_actual = EvaluateVector(solucion_actual,d);
        //Si la solución encontrada es mejor que la mejor solución hasta el momento, esta se convierte en la mejor
        if(evaluacion_actual < mejor_eval){
            mejor_solucion = solucion_actual;
            mejor_eval = evaluacion_actual;
        }
        else{
            //Si la mejor hasta el momento no es mejorada, esta pasa a ser la solucion actual y se reintenta la bl
            solucion_actual = mejor_solucion;
            evaluacion_actual = mejor_eval;
        }
    }
    return mejor_solucion;
}

void ES_Innicializado(Data &d, int max_evals, vector<int> &solucion_inicial,int &evaluaciones_iniciales){

    int evaluacion_entrada = EvaluateVector(solucion_inicial,d);
    int size = d.GetSize();

    int max_vecinos = 10*size;//Número máximo de vecinos que se pueden generar
    int max_exitos = floor(0.1*max_vecinos);//Número máximo de soluciones que se pueden aceptar
    vector<int> mejor_solucion; //Mejor solución obtenida hasta el momento
    vector<int> solucion_actual;//Solución de partida del problema y con la que se trabaja
    solucion_actual.resize(size);
    //GenerateRandomVector(solucion_actual);//Primera solución de la que partimos
    solucion_actual = solucion_inicial;
    //Mejor solución será igual a la primera genrada
    mejor_solucion = solucion_actual;
    int mejor_eval = EvaluateVector(mejor_solucion,d);
    bool exito = true; //Si en un enfriamiento(bucle interno) no se encuentra nada mejor, se pone a false
    int exitos = 0;
    int vecinos_generados = 0;
    int n_evals = evaluaciones_iniciales;
    int evaluacion_actual = mejor_eval;

    //Esquema de enfriamiento
    float mu = 0.4;
    float fi = 0.4;
    float t0 = mu*evaluacion_actual/-log(fi);
    float tf = pow(10,-3);
    int M = (float)max_evals/(float)max_vecinos;//Enfriamientos a realizar
    //cout << "M: " << M << endl;
    float B;//Variable necesaria para el enfriamiento de Cauchy
    B = (t0-tf)/(M*t0*tf);
    //cout << "t0: " << t0 << endl;

    float t = t0; //Temperatura del sistema
    int iteracion = 1; //Contador de iteraciones actuales
    int n_entradas_por_aleatorio = 0;


    while(n_evals < max_evals && exito){
        //Mientras el algoritmo mejore y el númer de evaluaciones sea menor al establecido

        exitos = 0;//Cuenta el númer de vecinos generados en cada enfriamiento
        for(int i = 0; i < max_vecinos && exitos < max_exitos;i++) {
            //Se genera de forma aleatoria el intercambio
            int p1;
            int p2;
            p1 = Randint(0,size-1);
            p2 = Randint(0,size-1);
            //para la pareja i-ésima
            int Swap = solucion_actual[p1];
            solucion_actual[p1] = solucion_actual[p2];
            solucion_actual[p2] = Swap;
            //Se comprueba la evaluación de la nueva solución creada
            int evaluacion_nueva = EvaluateNeigbour(evaluacion_actual, p1, p2, d, solucion_actual);
            n_evals++;
            //Ahora se decide si se acepta o si se rechaza la nueva solución
            if (evaluacion_nueva < evaluacion_actual) {
                //Si la nueva solución es mejor que la actual, se acepta sin más

                exitos++;

                //Ahora comprobamos si es la mejor solución encontrada
                if(evaluacion_nueva < mejor_eval){
                    mejor_solucion = solucion_actual;
                    mejor_eval = evaluacion_nueva;
                }
                evaluacion_actual = evaluacion_nueva;
            }
            else {
                //En este caso se aceptará una peor solución si la temperatura es lo bastante alta
                if (exp((-1 * (evaluacion_nueva - evaluacion_actual)) / (iteracion * t)) >= Randfloat(0, 1)) {

                    n_entradas_por_aleatorio++;

                    //Entonces puede escoger la nueva solución

                    exitos++;
                    if(evaluacion_nueva < mejor_eval){
                        mejor_solucion = solucion_actual;
                        mejor_eval = evaluacion_actual;
                    }
                    evaluacion_actual = evaluacion_nueva;
                }
                else{
                    //No escoge la solución actual
                    //Es necesario invertir los cambios realizados en el vector
                    Swap = solucion_actual[p1];
                    solucion_actual[p1] = solucion_actual[p2];
                    solucion_actual[p2] = Swap;
                }
            }
            vecinos_generados++;
        }
        if(exitos == 0){
            exito = false;
        }
        iteracion++;
        //Genera una nueva temperatura para el sistema
        t = t/(1+(B*t));
        //t = 0.99*t;

    }
    /*
    cout << "Numero de evaluaciones = " << n_evals << endl;
    cout << "numero de iteraciones = " << iteracion << endl;
    cout << "numero de entradas aleatorias = " << n_entradas_por_aleatorio << endl;
    */
    //evaluaciones_iniciales = n_evals;
    solucion_inicial = mejor_solucion;
}

vector<int>ILS_ES(Data &d, int max_evals, int n_iters){
    int size = d.GetSize();
    int size_sublista = size/4;

    vector<int> solucion_actual;
    solucion_actual.resize(size);
    GenerateRandomVector(solucion_actual);
    int evaluacion_actual = EvaluateVector(solucion_actual,d);
    int current_evals = 0;
    ES_Innicializado(d,max_evals,solucion_actual,current_evals);
    vector<int> mejor_solucion;
    mejor_solucion.resize(size);
    mejor_solucion = solucion_actual;
    int mejor_eval = evaluacion_actual;


    for(int i = 1; i < n_iters ;i++){
        //Se muta la solucion actual y se vuelve a lanzar la busqueda local.
        //Proceso de mutación
        //Primero se escoge el inicio de la sublista
        int inicio_sublista = Randint(0,size-size_sublista-1);

        //Se crea el vector de posiciones que serán alteradas a partir de la posición de inicio
        vector<int>posiciones;
        posiciones.resize(size_sublista);
        GenerateRandomVector(posiciones);
        for(int j = 0; j < size_sublista; j++){
            int swap = solucion_actual[j+inicio_sublista];
            solucion_actual[j+inicio_sublista] = solucion_actual[posiciones[j]];
            solucion_actual[posiciones[j]] = swap;
        }
        //BL sobre la solución mutada
        evaluacion_actual = EvaluateVector(solucion_actual,d);
        current_evals++;
        ES_Innicializado(d,max_evals,solucion_actual,current_evals);
        evaluacion_actual = EvaluateVector(solucion_actual,d);
        //Si la solución encontrada es mejor que la mejor solución hasta el momento, esta se convierte en la mejor
        if(evaluacion_actual < mejor_eval){
            mejor_solucion = solucion_actual;
            mejor_eval = evaluacion_actual;
        }
        else{
            //Si la mejor hasta el momento no es mejorada, esta pasa a ser la solucion actual y se reintenta la bl
            solucion_actual = mejor_solucion;
            evaluacion_actual = mejor_eval;
        }
    }
    return mejor_solucion;
}

vector<int>BMB(Data &d, int max_evals, int n_soluciones){
    int mejor_evaluacion;
    int current_eval;
    vector<int>solucion_actual;
    solucion_actual.resize(d.GetSize());
    GenerateRandomVector(solucion_actual);
    vector<int>mejor_solucion = solucion_actual;
    current_eval = EvaluateVector(solucion_actual,d);
    mejor_evaluacion = current_eval;
    int n_mejoras = 0;
    for(int i = 0; i < n_soluciones; i++){
        //Aplica una búsqueda local sobre la solución actual
        int iters = 0;
        LocalSearchInicializada(d,max_evals,solucion_actual,iters,current_eval);
        //Si la solución encontrada es la mejor, se actualiza la mejor solucion
        if(current_eval < mejor_evaluacion){
            mejor_evaluacion = current_eval;
            mejor_solucion = solucion_actual;
            n_mejoras++;
        }
        GenerateRandomVector(solucion_actual);
        current_eval = EvaluateVector(solucion_actual,d);
    }
    return mejor_solucion;
}


int main(int argc, char* argv[]) {
    if(argc != 3){
        cerr << "Se esperaba como argumento el nombre del archivo y la semilla\n";
        return 1;
    }
    Seed = atoi(argv[2]);
    cout << argv[1]<<" ";
    Data d;
    //d.LoadFromFile("/home/sergio/Escritorio/Informatica/3ºcurso/mh/practicas/p2/lipa90b.dat");
    d.LoadFromFile(argv[1]);
    vector<int> res;
    //Inicialización del cronómetro
    clock_t inicio;
    clock_t fin;


    //cout << "Enfriamiento simulado: ";
    inicio = clock();
    res = SimulatedAnnealing(d,50000);
    fin = clock();
    cout << EvaluateVector(res,d)<<" ";
    cout << ((float)fin - float(inicio))/CLOCKS_PER_SEC << endl;

    //cout << "Busqueda local: ";
    inicio = clock();
    res = LocalSearch(d,50000);
    fin = clock();
    cout << EvaluateVector(res,d)<<" ";
    cout << ((float)fin - float(inicio))/CLOCKS_PER_SEC << endl;

    //cout << "ILS: ";
    inicio = clock();
    res = ILS(d,50000,25);
    fin = clock();
    cout << EvaluateVector(res,d)<<" ";
    cout << ((float)fin - float(inicio))/CLOCKS_PER_SEC << endl;

    //cout << "ILS-ES: ";
    inicio = clock();
    res = ILS_ES(d,50000,25);
    fin = clock();
    cout << EvaluateVector(res,d)<<" ";
    cout << ((float)fin - float(inicio))/CLOCKS_PER_SEC << endl;

    //cout << "BMB: ";
    inicio = clock();
    res = BMB(d,50000,25);
    fin = clock();
    cout << EvaluateVector(res,d)<<" ";
    cout << ((float)fin - float(inicio))/CLOCKS_PER_SEC << endl;

    //cout << "GRASP: ";
    inicio = clock();
    res = Grasp(d,50000,25);
    fin = clock();
    cout << EvaluateVector(res,d) << " ";
    cout << ((float)fin - float(inicio))/CLOCKS_PER_SEC << endl;

    //cout << "Greedy: ";
    inicio = clock();
    res = greedy(d);
    fin = clock();
    cout << EvaluateVector(res,d) << " ";
    cout << ((float)fin - float(inicio))/CLOCKS_PER_SEC << endl;

}