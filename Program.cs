using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;
 



namespace Rextester
{

    public class Program
    {
        public static int[,] ReadFileCityToMatrix(string filename)
        {
            int[,] citymatrix;
            using (StreamReader sr = new StreamReader(filename))
            {
                string sfr = sr.ReadToEnd();
                string[] row = sfr.Split('\n');
                string[] col = row[0].Split(' ');
                citymatrix = new int[row.Length, col.Length];
                for(int i =0; i<row.Length; ++i)
                {
                    col = row[i].Split(' ');
                    for(int j = 0;j < col.Length; j++)
                    {
                        citymatrix[i, j] = Convert.ToInt32(col[j]);
                    }
                }
            }
            return citymatrix;
        }
        public static double[,] ReadFileCoordinatesCity(string filename)
        {
            double[,] citymatrix;
            using (StreamReader sr = new StreamReader(filename))
            {
                string sfr = sr.ReadToEnd();
                string[] row = sfr.Split('\n');
                string[] col = row[0].Split(' ');
                citymatrix = new double[row.Length, col.Length];
                for (int i = 0; i < row.Length; ++i)
                {
                    col = row[i].Split(' ');
                    for (int j = 0; j < col.Length; j++)
                    {
                        citymatrix[i, j] = Convert.ToInt32(col[j]);
                      
                    }
                  
                }
            }
            return citymatrix;
        }
        public static double Distance(int start, int end, double[,] coordinates)
        {
            double result = Math.Sqrt(Math.Pow((coordinates[end, 0] - coordinates[start, 0]),2) + Math.Pow((coordinates[end, 1] - coordinates[start, 1]) ,2));
            return result;
        }
        public static double Distance(double x1,double y1, double x2, double y2)
        {
            return Math.Sqrt(Math.Pow((x1 - x2), 2) + Math.Pow((y1 - y2), 2));
        }
        public static double[,] CoordinatesToDistanceMatrix(double[,] coordinates)
        {
            double[,] distanceMatrix = new double[coordinates.GetLength(0), coordinates.GetLength(0)];
            for(int i = 0; i < coordinates.GetLength(0); i++)
            {
                for(int j = 0; j < coordinates.GetLength(0); j++)
                {
                    if (i == j) distanceMatrix[i, j] = -1;
                    else
                    distanceMatrix[i, j] = Distance(i, j, coordinates);
                    Console.Write("{0:0.#} ", distanceMatrix[i, j]);
                }
                Console.WriteLine();
            }
            using (StreamWriter sw = new StreamWriter("242.txt", false, System.Text.Encoding.Default))
            {
                for (int i = 0; i < coordinates.GetLength(0); i++)
                {
                    for (int j = 0; j < coordinates.GetLength(0); j++)
                    {
                        sw.Write("{0:0.#} ", distanceMatrix[i,j]);
                    }
                    sw.WriteLine();
                }
            }
                return distanceMatrix;
        }
        public static void Main(string[] args)
        {
            //try
            //{
            //    Console.WriteLine("");
            string filename;
            Console.WriteLine("Enter filename: ");
            filename = Console.ReadLine();

            double[,] citycoordinates = ReadFileCoordinatesCity(filename);
            double[,] citymatrix = CoordinatesToDistanceMatrix(citycoordinates);
                //int [,] citymatrix= ReadFileCityToMatrix("1.txt");
            CitiesData citiesData = new CitiesData(citymatrix,citycoordinates);
            Console.WriteLine("Number of cities = " + citiesData.citiesDistance.GetLength(0));
            ClasterFogel cl = new ClasterFogel(5, citiesData);
            Stopwatch stopWatch = new Stopwatch();
            stopWatch.Start();

            cl.Clustering();
            Console.WriteLine("\n===================\nCLASTER ALGORITHM:");
            Console.WriteLine(cl.ToString());
            stopWatch.Stop();
            // Get the elapsed time as a TimeSpan value.
            TimeSpan ts = stopWatch.Elapsed;

            // Format and display the TimeSpan value.
            string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
                ts.Hours, ts.Minutes, ts.Seconds,
                ts.Milliseconds / 10);
            Console.WriteLine("RunTime " + elapsedTime);

            int totalNumberBees = 6;
                int numberInactive = 1;
                int numberActive = 3;
                int numberScout = 2;

                int maxNumberVisits = 3;
                int maxNumberCycles = 10;
           

                Hive hive = new Hive(totalNumberBees, numberInactive, numberActive, numberScout, maxNumberVisits, maxNumberCycles, citiesData);
               

                bool doProgressBar = false;
            Stopwatch stopWatch1 = new Stopwatch();
            stopWatch1.Start();
           
           
            hive.Solve(doProgressBar);

            Console.WriteLine("\n===================\nBee ALGORITHM:");
            Console.WriteLine(hive);
            stopWatch1.Stop();
            // Get the elapsed time as a TimeSpan value.
            TimeSpan ts1 = stopWatch1.Elapsed;

            // Format and display the TimeSpan value.
            string elapsedTime1 = String.Format("{0:00}:{1:00}:{2:00}.{3:00}",
                ts1.Hours, ts1.Minutes, ts1.Seconds,
                ts1.Milliseconds / 10);
            Console.WriteLine("RunTime " + elapsedTime1);
                Console.ReadKey();
             
        } // Main()
    } // class Program

    class Hive
    {
        public class Bee
        {
            public int status; // 0 = inactive, 1 = active, 2 = scout
            public int[] memoryMatrix; // problem-specific. a path of cities.
            public double measureOfQuality; // smaller values are better. total distance of path.
            public int numberOfVisits;

            public Bee(int status, int[] memoryMatrix, double measureOfQuality, int numberOfVisits)
            {
                this.status = status;
                this.memoryMatrix = new int[memoryMatrix.Length];
                Array.Copy(memoryMatrix, this.memoryMatrix, memoryMatrix.Length);
                this.measureOfQuality = measureOfQuality;
                this.numberOfVisits = numberOfVisits;
            }

            public override string ToString()
            {
                string s = "";
                s += "Status = " + this.status + "\n";
                s += " Memory = " + "\n";
                for (int i = 0; i < this.memoryMatrix.Length - 1; ++i)
                    s += this.memoryMatrix[i] + "->";
                s += this.memoryMatrix[this.memoryMatrix.Length - 1] + "\n";
                s += " Quality = " + this.measureOfQuality.ToString("F4");
                s += " Number visits = " + this.numberOfVisits;
                return s;
            }
        } // Bee

        static Random random = null; // multipurpose

        public CitiesData citiesData; // this is the problem-specific data we want to optimize

        public int totalNumberBees; // mostly for readability in the object constructor call
        public int numberInactive;
        public int numberActive;
        public int numberScout;

        public int maxNumberCycles; // one cycle represents an action by all bees in the hive
                                    //public int maxCyclesWithNoImprovement; // deprecated

        public int maxNumberVisits; // max number of times bee will visit a given food source without finding a better neighbor
        public double probPersuasion = 0.90; // probability inactive bee is persuaded by better waggle solution
        public double probMistake = 0.01; // probability an active bee will reject a better neighbor food source OR accept worse neighbor food source

        public Bee[] bees;
        public int[] bestMemoryMatrix; // problem-specific
        public double bestMeasureOfQuality;
        public int[] indexesOfInactiveBees; // contains indexes into the bees array

        public override string ToString()
        {
            string s = "";
            s += "Best path found: ";
            for (int i = 0; i < this.bestMemoryMatrix.Length - 1; ++i)
                s += this.bestMemoryMatrix[i] + "->";
            s += this.bestMemoryMatrix[this.bestMemoryMatrix.Length - 1] + "\n";

            s += "Path quality:    ";
            if (bestMeasureOfQuality < 10000.0)
                s += bestMeasureOfQuality.ToString("F4") + "\n";
            else
                s += bestMeasureOfQuality.ToString("#.####e+00");
            s += "\n";
            return s;
        }

        public Hive(int totalNumberBees, int numberInactive, int numberActive, int numberScout, int maxNumberVisits,
          int maxNumberCycles, CitiesData citiesData)
        {
            random = new Random(0);

            this.totalNumberBees = totalNumberBees;
            this.numberInactive = numberInactive;
            this.numberActive = numberActive;
            this.numberScout = numberScout;
            this.maxNumberVisits = maxNumberVisits;
            this.maxNumberCycles = maxNumberCycles;
            //this.maxCyclesWithNoImprovement = maxCyclesWithNoImprovement;

            //this.citiesData = new CitiesData(citiesData.cities.Length); // hive's copy of problem-specific data
            this.citiesData = citiesData; // reference to CityData

            // this.probPersuasion & this.probMistake are hard-coded in class definition

            this.bees = new Bee[totalNumberBees];
            this.bestMemoryMatrix =  citiesData.GenerateRandomMemoryMatrix(); // alternative initializations are possible
            this.bestMeasureOfQuality = citiesData.PathLength(this.bestMemoryMatrix);

            this.indexesOfInactiveBees = new int[numberInactive]; // indexes of bees which are currently inactive

            for (int i = 0; i < totalNumberBees; ++i) // initialize each bee, and best solution
            {
                int currStatus; // depends on i. need status before we can initialize Bee
                if (i < numberInactive)
                {
                    currStatus = 0; // inactive
                    indexesOfInactiveBees[i] = i; // curr bee is inactive
                }
                else if (i < numberInactive + numberScout)
                {
                    currStatus = 2; // scout
                }
                else
                {
                    currStatus = 1; // active
                }

                int[] randomMemoryMatrix =citiesData.GenerateRandomMemoryMatrix();
                double mq = citiesData.PathLength(randomMemoryMatrix);
                int numberOfVisits = 0;

                bees[i] = new Bee(currStatus, randomMemoryMatrix, mq, numberOfVisits); // instantiate current bee
                
                // does this bee have best solution?
                if (bees[i].measureOfQuality < bestMeasureOfQuality) // curr bee is better (< because smaller is better)
                {
                    Array.Copy(bees[i].memoryMatrix, this.bestMemoryMatrix, bees[i].memoryMatrix.Length);
                     this.bestMeasureOfQuality = bees[i].measureOfQuality;
                }
                
            } // each bee

        } // TravelingSalesmanHive ctor






        public void Solve(bool doProgressBar) // find best Traveling Salesman Problem solution
        {
            bool pb = doProgressBar; // just want a shorter variable
            int numberOfSymbolsToPrint = 10; // 10 units so each symbol is 10.0% progress
            int increment = this.maxNumberCycles / numberOfSymbolsToPrint;
            if (pb) Console.WriteLine("\nEntering SBC Traveling Salesman Problem algorithm main processing loop\n");
            if (pb) Console.WriteLine("Progress: |==========|"); // 10 units so each symbol is 10% progress
            if (pb) Console.Write("           ");
            int cycle = 0;
            using (StreamWriter sw = new StreamWriter("232.txt", false, System.Text.Encoding.Default))
            {
                while (cycle < this.maxNumberCycles)
                {
                    for (int i = 0; i < this.totalNumberBees; ++i) // each bee
                    {

                        sw.WriteLine(bees[i].ToString());


                        if (this.bees[i].status == 1) // active bee
                            ProcessActiveBee(i);
                        else if (this.bees[i].status == 2) // scout bee
                            ProcessScoutBee(i);
                        else if (this.bees[i].status == 0) // inactive bee
                            ProcessInactiveBee(i);

                    } // for each bee
                    ++cycle;
                    sw.WriteLine();
                    // print a progress bar
                    if (pb && cycle % increment == 0)
                        Console.Write("^");
                } // main while processing loop

                if (pb) Console.WriteLine(""); // end the progress bar
            }
        }// Solve()

        private void ProcessInactiveBee(int i)
        {
            return; // not used in this implementation
        }

        private void ProcessActiveBee(int i)
        {
            int[] neighbor = citiesData.GenerateNeighborMemoryMatrix(bees[i].memoryMatrix); // find a neighbor solution
            double neighborQuality = citiesData.PathLength(neighbor); // get its quality
            double prob = random.NextDouble(); // used to determine if bee makes a mistake; compare against probMistake which has some small value (~0.01)
            bool memoryWasUpdated = false; // used to determine if bee should perform a waggle dance when done
            bool numberOfVisitsOverLimit = false; // used to determine if bee will convert to inactive status

            if (neighborQuality < bees[i].measureOfQuality) // active bee found better neighbor (< because smaller values are better)
            {
                if (prob < probMistake) // bee makes mistake: rejects a better neighbor food source
                {
                    ++bees[i].numberOfVisits; // no change to memory but update number of visits
                    if (bees[i].numberOfVisits > maxNumberVisits) numberOfVisitsOverLimit = true;
                }
                else // bee does not make a mistake: accepts a better neighbor
                {
                    Array.Copy(neighbor, bees[i].memoryMatrix, neighbor.Length); // copy neighbor location into bee's memory
                    bees[i].measureOfQuality = neighborQuality; // update the quality
                    bees[i].numberOfVisits = 0; // reset counter
                    memoryWasUpdated = true; // so that this bee will do a waggle dance 
                }
            }
            else // active bee did not find a better neighbor
            {
                //Console.WriteLine("c");
                if (prob < probMistake) // bee makes mistake: accepts a worse neighbor food source
                {
                    Array.Copy(neighbor, bees[i].memoryMatrix, neighbor.Length); // copy neighbor location into bee's memory
                    bees[i].measureOfQuality = neighborQuality; // update the quality
                    bees[i].numberOfVisits = 0; // reset
                    memoryWasUpdated = true; // so that this bee will do a waggle dance 
                }
                else // no mistake: bee rejects worse food source
                {
                    ++bees[i].numberOfVisits;
                    if (bees[i].numberOfVisits > maxNumberVisits) numberOfVisitsOverLimit = true;
                }
            }

            if (numberOfVisitsOverLimit == true)
            {
                bees[i].status = 0; // current active bee transitions to inactive
                bees[i].numberOfVisits = 0; // reset visits (and no change to this bees memory)
                int x = random.Next(numberInactive); // pick a random inactive bee. x is an index into a list, not a bee ID
                bees[indexesOfInactiveBees[x]].status = 1; // make it active
                indexesOfInactiveBees[x] = i; // record now-inactive bee 'i' in the inactive list
            }
            else if (memoryWasUpdated == true) // current bee returns and performs waggle dance
            {
                // first, determine if the new memory is a global best. note that if bee has accepted a worse food source this can't be true
                if (bees[i].measureOfQuality < this.bestMeasureOfQuality) // the modified bee's memory is a new global best (< because smaller is better)
                {
                    Array.Copy(bees[i].memoryMatrix, this.bestMemoryMatrix, bees[i].memoryMatrix.Length); // update global best memory
                    this.bestMeasureOfQuality = bees[i].measureOfQuality; // update global best quality
                }
                DoWaggleDance(i);
            }
            else // number visits is not over limit and memory was not updated so do nothing (return to hive but do not waggle)
            {
                return;
            }
        } // ProcessActiveBee()

        private void ProcessScoutBee(int i)
        {
            int[] randomFoodSource = citiesData.GenerateRandomMemoryMatrix(); // scout bee finds a random food source. . . 
            double  randomFoodSourceQuality = citiesData.PathLength(randomFoodSource); // and examines its quality
            if (randomFoodSourceQuality < bees[i].measureOfQuality) // scout bee has found a better solution than its current one (< because smaller measure is better)
            {
                Array.Copy(randomFoodSource, bees[i].memoryMatrix, randomFoodSource.Length); // unlike active bees, scout bees do not make mistakes
                bees[i].measureOfQuality = randomFoodSourceQuality;
                // no change to scout bee's numberOfVisits or status

                // did this scout bee find a better overall/global solution?
                if (bees[i].measureOfQuality < bestMeasureOfQuality) // yes, better overall solution (< because smaller is better)
                {
                    Array.Copy(bees[i].memoryMatrix, this.bestMemoryMatrix, bees[i].memoryMatrix.Length); // copy scout bee's memory to global best
                    this.bestMeasureOfQuality = bees[i].measureOfQuality;
                } // better overall solution

                DoWaggleDance(i); // scout returns to hive and does waggle dance

            } // if scout bee found better solution
        } // ProcessScoutBee()

        private void DoWaggleDance(int i)
        {
            for (int ii = 0; ii < numberInactive; ++ii) // each inactive/watcher bee
            {
                int b = indexesOfInactiveBees[ii]; // index of an inactive bee
                if (bees[b].status != 0) throw new Exception("Catastrophic logic error when scout bee waggles dances");
                if (bees[b].numberOfVisits != 0) throw new Exception("Found an inactive bee with numberOfVisits != 0 in Scout bee waggle dance routine");
                if (bees[i].measureOfQuality < bees[b].measureOfQuality) // scout bee has a better solution than current inactive/watcher bee (< because smaller is better)
                {
                    double p = random.NextDouble(); // will current inactive bee be persuaded by scout's waggle dance?
                    if (this.probPersuasion > p) // this inactive bee is persuaded by the scout (usually because probPersuasion is large, ~0.90)
                    {
                        Array.Copy(bees[i].memoryMatrix, bees[b].memoryMatrix, bees[i].memoryMatrix.Length);
                        bees[b].measureOfQuality = bees[i].measureOfQuality;
                    } // inactive bee has been persuaded
                } // scout bee has better solution than watcher/inactive bee
            } // each inactive bee
        } // DoWaggleDance()

    } // class ShortestPathHive

    class CitiesData
    {
        bool mod;
        int start, end;
        public int[] cities;
        public double[,] citiesDistance;
        public double[,] coordinates;
        //public int start;
        private Random rand = new Random(0);
        public CitiesData(double[,] citiesDistance,double[,]coordinates,bool mod=false,int start = 0, int end = 0,int[]cities=null)
        {
            this.cities = cities;
            this.start = start;
            this.end = end;
            this.mod = mod;
            this.coordinates = coordinates;
            this.citiesDistance = citiesDistance; 
           // start = startCity;
        }
        public double PathLength(int[] path)
        {
            double answer = 0;
            for (int i = 0; i < path.Length-1; ++i)
            {
               
                answer += citiesDistance[path[i],path[i+1]];
            }
            answer += citiesDistance[path[path.Length-1], path[0]];                                                                                                               
            return answer;
        }
        public int[] GenerateRandomMemoryMatrix()
        {
            int[] result;
             // // problem-specific
            if (!mod)
            {
                result = new int[this.citiesDistance.GetLength(0)];
                result[0] = 0;
                for (int i = 1; i < this.citiesDistance.GetLength(0); ++i)
                {
                    List<int> arage = new List<int>();
                    for (int j = 0; j < this.citiesDistance.GetLength(0); ++j)
                    {
                        if (citiesDistance[result[i - 1], j] != -1)
                        {
                            bool flag = true;
                            for (int k = 0; k < i; ++k)
                            {
                                if (result[k] == j) flag = false;
                            }
                            if (flag) arage.Add(j);
                        }
                    }
                    if (arage.Count == 0) throw new Exception("no path");
                    int randomcity = rand.Next(0, arage.Count);
                    result[i] = arage[randomcity];
                }
            }
            else {
                bool f = false;
                 result = new int[this.cities.Length];
                result[0] = start;
                if (start != end) f = true;
                if(f)result[this.cities.Length-1] = end;
                int t = 0;
                if (!f) {
                    for (int i = 1; i < cities.Length; ++i)
                    {
                        List<int> arage = new List<int>();
                        foreach (int j in cities)
                        {
                            if (citiesDistance[result[i - 1], j] != -1)
                            {
                                bool flag = true;
                                for (int k = 0; k < i; ++k)
                                {
                                    if (result[k] == j) flag = false;
                                }
                                if (flag) arage.Add(j);
                            }
                        }
                        if (arage.Count == 0) throw new Exception("no path");
                        int randomcity = rand.Next(0, arage.Count);
                        result[i] = arage[randomcity];
                    }
                }
                else
                {
                    for (int i = 1; i < cities.Length - 1; ++i)
                    {
                        List<int> arage = new List<int>();
                        foreach (int j in cities)
                        {
                            if (citiesDistance[result[i - 1], j] != -1)
                            {
                                bool flag = true;
                                for (int k = 0; k < i; ++k)
                                {
                                    if (result[k] == j) flag = false;
                                }
                                if (flag) arage.Add(j);
                            }
                        }
                        if (arage.Count == 0) throw new Exception("no path");
                        int randomcity = rand.Next(0, arage.Count);
                        result[i] = arage[randomcity];
                    }
                }
            }
            return result;
        } // GenerateRandomMemoryMatrix()

        public int[] GenerateNeighborMemoryMatrix(int[] memoryMatrix)
        {
            if (memoryMatrix.Length == 1) return memoryMatrix;
            int[] result = new int[memoryMatrix.Length];
            Array.Copy(memoryMatrix, result, memoryMatrix.Length);

            int ranIndex = rand.Next(1, result.Length-1); // [1, Length-2] inclusive
           
            int adjIndex;
            if (ranIndex == result.Length - 1)
                adjIndex = 1;
            else      
                adjIndex = ranIndex + 1;

            int tmp = result[ranIndex];
            result[ranIndex] = result[adjIndex];
            result[adjIndex] = tmp;

            return result;
        } // GenerateNeighborMemoryMatrix()
        
    } // class CitiesData



    class ClasterFogel
    {
        private int R;
        public int[] resultPath;
        public List<Cluster> clasters;
        private CitiesData citiesData;
        public ClasterFogel(int R,CitiesData citiedata)
        {
            this.citiesData = citiedata;
            this.R = R;
        }
        public class Cluster
        {
            public int inputPathTo;
            public int inputPathFrom;
            public int outputPathTo;
            public int outputPathFrom;
            public int inputClaster;
            public int outputClaster;
            public int[] cities;
            public int[] path;
            public double[] center;
            public void setInput(int it, int inpf)
            {
                inputPathTo = it;
                inputPathFrom = inpf;

            }
            public void setOutput(int ot, int of)
            {
                outputPathFrom = of;
                outputPathTo = ot;

            }
            public void setClusterInput(int inp)
            {
                inputClaster = inp;
            }
            public void setClusterOutp(int outp)
            {
                outputClaster = outp;
            }
        }
       
        private double[] CenterMass(List<int> cities)
        {
            double[] xy = new double[2];
            xy[0] = xy[1] = 0;
            //double P = 0;
            for (int i = 0; i <cities.Count; i++)
            {
                // применяем формулы (*)
                
                xy[0] +=  citiesData.coordinates[cities[i],0];
                xy[1] += citiesData.coordinates[cities[i], 1];


            }
            xy[0] /= cities.Count;
            xy[1] /= cities.Count;
            return xy;
        }
        public void Clustering()
        {
            Random rand = new Random();
            //bool[] Pool = new bool[citiedata.citiesDistance.GetLength(0)];
            List<int> poolCities = new List<int>();
            for (int i = 0; i < citiesData.citiesDistance.GetLength(0); i++)
            {
                poolCities.Add(i);
            }

            List<Cluster> Clusters = new List<Cluster>();
            do
            {

                double[] centerCluster = new double[2];
                int city = poolCities[rand.Next(poolCities.Count)];
                centerCluster[0]=citiesData.coordinates[city,0];
                centerCluster[1] = citiesData.coordinates[city, 1];
                double[] newcenterCluster = new double[2] ;
                List<int> citiesFromCluster = new List<int>();
                do
                {
                    citiesFromCluster.Clear();
                    newcenterCluster = centerCluster;
                    foreach(int i in poolCities)
                    {
                           if(Program.Distance(centerCluster[0],centerCluster[1],citiesData.coordinates[i,0], citiesData.coordinates[i, 1]) <= R)
                        {
                            citiesFromCluster.Add(i);
                        }
                    }
                    centerCluster = CenterMass(citiesFromCluster);
                } while (!(newcenterCluster[0]==centerCluster[0]&&newcenterCluster[1]==centerCluster[1]));
                foreach(int i in citiesFromCluster)
                {
                    poolCities.Remove(i);
                }
                Clusters.Add(new Cluster() { cities = citiesFromCluster.ToArray(), center = newcenterCluster });

            } while (poolCities.Count != 0);
            //Clustering
            //Start SearchWay
            int totalNumberBees = 6;
            int numberInactive = 1;
            int numberActive = 3;
            int numberScout = 2;

            int maxNumberVisits = 3;
            int maxNumberCycles = 10;


            double[,] ncitiesCenters = new double[Clusters.Count,2];
            for(int i=0;i<Clusters.Count;i++)
            {
                ncitiesCenters[i,0] = Clusters[i].center[0];
                ncitiesCenters[i, 1] = Clusters[i].center[1];
            }
            CitiesData ncd = new CitiesData(Program.CoordinatesToDistanceMatrix(ncitiesCenters), ncitiesCenters);
            Hive hive = new Hive(totalNumberBees, numberInactive, numberActive, numberScout, maxNumberVisits, maxNumberCycles, ncd);
          

            bool doProgressBar = false;
            hive.Solve(doProgressBar);
           
            for(int i = 0; i < Clusters.Count-1; i++)
            {
                ShortPathClusters(Clusters[hive.bestMemoryMatrix[i]], Clusters[hive.bestMemoryMatrix[i + 1]]);
                
                
            }
            ShortPathClusters(Clusters[Clusters.Count-1], Clusters[0]);
          
            foreach (Cluster c in Clusters)
            {
                ShortPathInsideCluster(c);
            }
            int[] result = ResultPath(Clusters,hive);
            resultPath = result;
            this.clasters = Clusters;
            
        }
        public override string ToString()
        {
            string res = "";
            res += "Clasters: ";
            foreach(Cluster c in clasters)
            {
                res += "\n";
                foreach(int i in c.cities)
                {
                    res += i;
                }
            }
            res += "\n Path: ";
            foreach(int i in resultPath)
            {
                res += "->"+i;
            }
            res += "\n Quality:" + citiesData.PathLength(resultPath);
            return res;
        }
     
        private int[]ResultPath(List<Cluster> c,Hive hive)
        {
            int[] path = new int[citiesData.citiesDistance.GetLength(0)];
            int t = 0;
            bool[] clastersFlag = new bool[c.Count];
            for(int i = 0; i < clastersFlag.Length; i++)
            {
                clastersFlag[i] = false;
            }
            foreach(int i in hive.bestMemoryMatrix)
            {
               
                for(int j=0;j<c[i].path.Length-1;j++)
                {
                    
                    path[t] = c[i].path[j];
                    t++;
                  
                }
                if (c[i].path[c[i].path.Length - 1] != c[i].outputPathFrom)
                {
                    path[t] = c[i].path[c[i].path.Length-1];
                    t++;
                }


            }
            return path;

            
        }
        private void ShortPathInsideCluster(Cluster c)
        {
            int totalNumberBees = 6;
            int numberInactive = 1;
            int numberActive = 3;
            int numberScout = 2;

            int maxNumberVisits = 3;
            int maxNumberCycles = 10;



            CitiesData ncdc = new CitiesData(citiesData.citiesDistance, citiesData.coordinates, true, c.inputPathTo, c.outputPathFrom, c.cities);
            Hive hive = new Hive(totalNumberBees, numberInactive, numberActive, numberScout, maxNumberVisits, maxNumberCycles, ncdc);
            hive.Solve(false);
            c.path = hive.bestMemoryMatrix;

        }
  

        private void ShortPathClusters(Cluster c1, Cluster c2)
        {
            double shrt = 999999999;
            int ct1=-1;
            int ct2=-1;
            foreach(int city in c1.cities)
            {
                foreach(int city1 in c2.cities)
                {
                    if (citiesData.citiesDistance[city, city1] < shrt)
                    {
                        shrt = citiesData.citiesDistance[city, city1];
                        ct1 = city;
                        ct2 = city1;

                    }
                }
            }
            c1.setOutput(ct2, ct1);
            c2.setInput(ct2, ct1);
           

            
        }

    }

} 

