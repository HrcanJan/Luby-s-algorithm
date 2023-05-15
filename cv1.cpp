#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <string>
#include <algorithm>
#include <time.h>
#include <fstream>
#include <random>
#include <mpi.h>
#include <chrono>
#include <iomanip>
#include <fstream>

using namespace std;

float random_num(int d, int rank) {
    int random;
    srand(time(NULL));

    if (d >= 1) {
        random = rand() % static_cast<int>(pow(2, d)) + 1;
        cout << rank << ": " << random << endl;
        if (random == 1)
            return 1;
        return 0;
    }
    return 1;
}


vector<vector<int>> greedyPartition(vector<vector<int>> numbers, int partitions, int n) {

    int num1 = ceil(static_cast<double>(n) / partitions);
    int i = 0;
    vector<int> vertices(n);

    for (int i = 0; i < n; i++)
        vertices[i] = i;

    vector<vector<int>> result(partitions);

    while (vertices.size() > 0) {

        if (vertices.size() > num1) {

            int u = vertices[0];
            while ((result[i].size() + numbers[u].size() + 1) <= num1) {

                result[i].push_back(u);
                for (int k = 0; k < numbers[u].size(); k++) {

                    auto iter = find(vertices.begin(), vertices.end(), numbers[u][k]);
                    if (iter != vertices.end()) {
                        result[i].push_back(numbers[u][k]);
                    }
                }

                auto iter = find(vertices.begin(), vertices.end(), u);
                if (iter != vertices.end()) {
                    vertices.erase(iter);
                }

                for (int k = 0; k < numbers[u].size(); k++) {
                    auto iter = find(vertices.begin(), vertices.end(), numbers[u][k]);
                    if (iter != vertices.end()) {
                        vertices.erase(iter);
                    }
                }

                u = vertices[0];

            }
            while ((result[i].size() + 1) <= num1) {

                result[i].push_back(u);
                auto iter = find(vertices.begin(), vertices.end(), u);
                if (iter != vertices.end()) {
                    vertices.erase(iter);
                }
                u = vertices[0];
            }

            ++i;
        }
        else {

            for (int k = 0; k < vertices.size(); k++)
                result[i].push_back(vertices[k]);

            vertices.clear();
        }
    }
    return result;
}


vector<int> luby(const vector<vector<int>>& adjList, vector<int> nodes) {
    // Initialize the set of unmarked nodes
    vector<int> unmarked(nodes.size());
    unmarked = nodes;

    // Initialize the set of marked nodes
    vector<int> marked(nodes.size());

    // Initialize the set of independent nodes
    vector<int> independent;

    // Repeat until all nodes are marked
    while (!unmarked.empty()) {

        // Choose a random node to mark
        int x = rand() % unmarked.size();
        int u = unmarked[x];
        auto y = find(nodes.begin(), nodes.end(), u);
        int z = distance(nodes.begin(), y);
        marked[z] = 1;
        unmarked.erase(find(unmarked.begin(), unmarked.end(), u));

        // Check if u is independent
        bool is_independent = true;
        for (int i = 0; i < adjList[u].size(); i++) {
            int v = adjList[u][i];

            auto it = find(nodes.begin(), nodes.end(), v);
            if (it != nodes.end()) {
                int index = distance(nodes.begin(), it);
                if (index < marked.size() && marked[index]) {
                    is_independent = false;
                    break;
                }
            }
        }

        // If u is independent, add it to the independent set
        if (is_independent) {
            independent.push_back(u);

            // Remove u's neighbors from the unmarked set
            vector<int> unmarked_copy = unmarked;
            for (int i = 0; i < adjList[u].size(); i++) {
                int v = adjList[u][i];
                if (find(unmarked_copy.begin(), unmarked_copy.end(), v) != unmarked_copy.end()) {
                    unmarked.erase(remove(unmarked.begin(), unmarked.end(), v), unmarked.end());
                }
            }
        }

    }

    return independent;
}

int main(int argc, char** argv) {

    // Start the timer
    auto start = std::chrono::high_resolution_clock::now();

    // Initialize MPI
    MPI_Init(&argc, &argv);
    int rank, num;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num);

    srand(time(NULL));

    // Open file for reading

    vector<int> nodes;
    vector<vector<int>> adjList;
    ifstream inputFile("data.txt");

    if (inputFile.is_open()) {

        int u, v;
        while (inputFile >> u >> v) {

            if (u >= adjList.size()) {
                adjList.resize(u + 1);
            }

            if (v >= adjList.size()) {
                adjList.resize(v + 1);
            }

            adjList[u].push_back(v);
            adjList[v].push_back(u);
            nodes.push_back(u);
            nodes.push_back(v);
        }

        inputFile.close();
    }
    else {

        cout << "Unable to open input file" << endl;
    }

    int n = adjList.size();
    sort(nodes.begin(), nodes.end());                       // Sort the vector
    auto last = std::unique(nodes.begin(), nodes.end());    // Remove consecutive duplicates
    nodes.erase(last, nodes.end());                         // Erase redundant elements

    // Initialize the sparse graph
    // int n = 11;  // Number of nodes
    // int m = 11;  // Number of edges
    // vector<int> nodes = { 0,1,2,3,4,5,6,7,8,9,10 };
    // vector<vector<int>> adjList(n);
    // for (int i = 0; i < m; i++) {
    //    int u = rand() % n;
    //    int v = rand() % n;
    //    adjList[u].push_back(v);
    //    adjList[v].push_back(u);
    //}

    //adjList[0].push_back(1);
    //adjList[1].push_back(0);
    //adjList[0].push_back(2);
    //adjList[2].push_back(0);
    //adjList[1].push_back(2);
    //adjList[2].push_back(1);
    //adjList[1].push_back(3);
    //adjList[3].push_back(1);
    //adjList[3].push_back(4);
    //adjList[4].push_back(3);
    //adjList[2].push_back(4);
    //adjList[4].push_back(2);
    //adjList[5].push_back(1);
    //adjList[1].push_back(5);
    //adjList[5].push_back(6);
    //adjList[6].push_back(5);
    //adjList[6].push_back(4);
    //adjList[4].push_back(6);
    //adjList[7].push_back(4);
    //adjList[4].push_back(7);
    //adjList[6].push_back(7);
    //adjList[7].push_back(6);

    double aa = log10(n);
    int partitions = ceil(aa);

    if (rank == 0) {
        // Print the adjacency list
        //for (int i = 0; i < adjList.size(); i++) {
        //    cout << "Node " << i << ": ";
        //    for (int j = 0; j < adjList[i].size(); j++) {
        //        cout << adjList[i][j] << " ";
        //    }
        //    cout << endl;
        //}

        //for (int i = 0; i < nodes.size(); i++) {
        //    cout << "Node " << nodes[i] << endl;
        //}
        cout << n << " nodes." << endl << partitions << " partitions." << endl;
    }

    
    vector<vector<int>> partitionsResult = greedyPartition(adjList, partitions, n);
    vector<int> finale;

    // Distribute the work among MPI processes
    int startIdx = rank * (partitionsResult.size() / num);
    int endIdx = (rank + 1) * (partitionsResult.size() / num);


    MPI_Barrier(MPI_COMM_WORLD);
    // Print the partitions
    for (int i = startIdx; i < endIdx; i++) {

        // cout << "Partition " << i + 1 << ": ";
        vector<vector<int>> pom;
        vector<int> pom_nodes;

        for (int j = 0; j < partitionsResult[i].size(); j++) {
            // cout << partitionsResult[i][j] << " ";
            pom_nodes.push_back(partitionsResult[i][j]);
            pom.push_back(adjList[partitionsResult[i][j]]);
        }
        // cout << endl;

        vector<int> test = luby(adjList, pom_nodes);
        // cout << "Independent set: ";
        for (int k = 0; k < test.size(); k++) {
            // cout << test[k] << " ";
            finale.push_back(test[k]);
        }
        // cout << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Gather the sizes of the finale vectors from all processes
    vector<int> sizes(num);
    int finaleSize = finale.size();
    MPI_Gather(&finaleSize, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Calculate displacement array for MPI_Gatherv
    vector<int> displacements(num);
    displacements[0] = 0;
    for (int i = 1; i < num; i++) {
        displacements[i] = displacements[i - 1] + sizes[i - 1];
    }

    // Gather the finale vectors from all processes
    vector<int> gatheredFinale(displacements[num - 1] + sizes[num - 1]);
    MPI_Gatherv(finale.data(), finaleSize, MPI_INT, gatheredFinale.data(), sizes.data(), displacements.data(), MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        vector<int> independent = luby(adjList, gatheredFinale);

        // Print the independent set
        cout << "Independent set: ";
        for (int i = 0; i < independent.size(); i++) {
            cout << independent[i] << " ";
        }
        cout << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Stop the timer
    auto end = chrono::high_resolution_clock::now();

    // Calculate the elapsed time
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);

    // Convert duration to seconds
    double seconds = duration.count() / 1000.0;
    vector<double> allTimes(num);
    MPI_Gather(&seconds, 1, MPI_DOUBLE, allTimes.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // Process 0 can now access all the times in the 'allTimes' vector
        for (int i = 0; i < num; i++) {
            cout << "Compilation time from process " << i << ": " << allTimes[i] << " seconds" << endl;
        }
    }

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
