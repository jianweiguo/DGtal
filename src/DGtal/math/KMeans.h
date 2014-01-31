#include <iostream>
#include <vector>
#include <stdlib.h>


template<typename Point, typename DistFunctor>
void KMeans(const std::vector<Point> &S,
	    const unsigned int k,
	    const DistFunctor &D,
	    std::vector<unsigned int> &registration,
	    std::vector<Point> &Centroid)
{
  //Init
  std::vector<Point> Sp = S;
  Centroid.resize(k); //Cluster centroids
  registration.resize(S.size()); //For each point, the cluster Id

  //We create a first set of centroid
  for(unsigned int i=0; i < k; i++)
    {
      unsigned int j = rand() % Sp.size();
      bool alreadySet = false;
      do
        {
          alreadySet = false;
          for( unsigned int ii = 0; ii < i; ++ii )
            {
              if( Centroid[ii] == Sp[j])
                {
                  alreadySet = true;
                  break;
                }
            }
          if( alreadySet )
            {
              j = rand() % Sp.size();
            }
        }
      while( alreadySet );
      Centroid[i] = Sp[j];
      Sp[j]  = Sp[Sp.size() - 1];
      Sp.pop_back();
    }
  Sp.clear();

  std::cout<<"   first assign"<<std::endl;

  //First assignment
  for(unsigned int i = 0; i < S.size(); i++)
    {
      registration[i] = 0;
      for(unsigned int c =  1 ; c < k; c++)
        if (D.distance(S[i], Centroid[c])< D.distance(S[i], Centroid[registration[i]]))
          registration[i] = c;
    }
  
  std::cout<<"   first centroids"<<std::endl;
  //Centroids
  for(unsigned int i = 0 ; i < k; ++i)
    {
      Point meanPoint= Centroid[i];
      unsigned int count=1;
      for(unsigned int j = 0; j < S.size(); j++)
        if (registration[j] == i)
          {
            meanPoint += S[j];
            count++;
          }
      Centroid[i] = meanPoint/(double)count;
    }

  
  //Main loop
  bool change = true;
  
  std::cout<<" Main loop: "<<std::endl;
  while (change)
    {
      std::cout << "#";
      change = false;
      unsigned int id;
      //Reassignment
      for(unsigned int i = 0; i < S.size(); i++)
        {
          id = registration[i];
          //	  std::cout << "Si="<<S[i].transpose()<<"  registration: "<<id
          //		    << "  min="<<Centroid[id].transpose()
          //		    << " d="<<D.distance(S[i], Centroid[id])<<std::endl;
          for(unsigned int c = 0 ; c < k; c++)
            {
              //  std::cout << " dist "<<S[i].transpose()<< " -- "<<Centroid[c].transpose()
              //<< "   d = "<< D.distance(S[i], Centroid[c])<<std::endl;
              if (D.distance(S[i], Centroid[c]) < D.distance(S[i], Centroid[id]))
                id = c;
            }
          //	  std::cout << "   id="<<id<<std::endl;
          if (id != registration[i])
            {
              //update detected
              registration[i] = id;
              change = true;
            }

        }
      
      //Centroid relaxation
      if (change)
        {
          for(unsigned int i = 0 ; i < k; ++i)
            {
              Point meanPoint = Centroid[i];
              unsigned int count=1;
              for(unsigned int j = 0; j < S.size(); j++)
                if (registration[j] == i)
                  {
                    meanPoint += S[j];
                    count ++;
                  }
              Centroid[i] = meanPoint/(double)count;
            }
        }

      //    /  for(unsigned int i=0; i <k ; i++)
      //	std::cout << Centroid[i].transpose()<<std::endl;

      
    }
  std::cout<<std::endl;
}
