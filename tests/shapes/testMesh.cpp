/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

/**
 * @file test3DMeshFromPoints.cpp
 * @ingroup Tests
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2012/07/02
 *
 * Functions for testing class 3DMeshFromPoints.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/writers/MeshWriter.h"
#include "DGtal/shapes/Mesh.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace DGtal::Z2i;



///////////////////////////////////////////////////////////////////////////////
// Functions for testing class Mesh.
///////////////////////////////////////////////////////////////////////////////

/**
 * Test the mesh object construction.
 */
bool testMesh()
{
  
  trace.beginBlock ( "Testing Mesh  ..." );
  bool ok = true;
  trace.beginBlock ( "Testing Mesh contruction  ..." );
  Mesh<Point> aMesh;
  Point p0=Point(0,0);
  Point p1=Point(0,1);
  Point p2=Point(1,2);
  Point p3=Point(3,2);
  Point p4=Point(3,3);
  Point p5=Point(3,4);
  Point p6=Point(4,6);
  
  aMesh.addVertex(p0);
  aMesh.addVertex(p1);
  aMesh.addVertex(p2);
  aMesh.addVertex(p3);
  aMesh.addVertex(p4);
  aMesh.addVertex(p5);
  
  aMesh.addTriangularFace(0,1,2);
  aMesh.addTriangularFace(3,4,5);
  
  Mesh<Point>::MeshFace tface0 = aMesh.getFace(0);
  Mesh<Point>::MeshFace tface1 = aMesh.getFace(1);
  Point p0f0 = aMesh.getVertex(tface0.at(0));
  Point p1f0 = aMesh.getVertex(tface0.at(1));
  Point p2f0 = aMesh.getVertex(tface0.at(2));

  Point p0f1 = aMesh.getVertex(tface1.at(0));
  Point p1f1 = aMesh.getVertex(tface1.at(1));
  Point p2f1 = aMesh.getVertex(tface1.at(2));
  trace.info() << "Set of points" << endl;
  trace.info() << p0 << p1 << p2 << endl;
  trace.info() << p3 << p4 << p5 << endl;
  
  trace.info() << "Face1 points " << endl;
  trace.info() << p0f0 << p1f0 << p2f0<< endl;

  trace.info() << "Face2 points " << endl;
  trace.info() << p0f1 << p1f1 << p2f1<< endl;

  
  bool okMeshConstruct =  (p0==p0f0) && (p1==p1f0) && (p2==p2f0) && 
    (p3==p0f1) && (p4==p1f1) && (p5==p2f1) ;
  
  
  trace.endBlock();
  bool okMeshIterators = true;
  trace.beginBlock ( "Testing Mesh iterator  ..." );
  unsigned int nb=0;
  // just testing nb iterations on const iterator
  for(  Mesh<Point>::VertexStorage::const_iterator it = aMesh.cVertexBegin(); 
       it !=aMesh.cVertexEnd(); 
       it++){
    nb++;    
  }
  okMeshIterators = nb == aMesh.nbVertex();
  // testing to change vertex  on  iterator
  for(  Mesh<Point>::VertexStorage::iterator it = aMesh.vertexBegin(); 
       it !=aMesh.vertexEnd(); 
       it++){
    (*it)[0]+=10.0; (*it)[1]+=5.0;
  }
  okMeshIterators = okMeshIterators &&  (aMesh.getVertex(5))[0]==13;

  // testing changing color of individual face:
  aMesh.setFaceColor(1, DGtal::Color::Red);
  bool okMeshColor = (aMesh.getFaceColor(0)==DGtal::Color::White)
                     && (aMesh.getFaceColor(1)==DGtal::Color::Red) ;
  
  ok = ok & okMeshConstruct &&  okMeshIterators && okMeshColor;  
  trace.endBlock();
  return ok;

}



/**
 * Test mesh generation from static fonctions.
 */
bool testMeshGeneration()
{
  
  trace.beginBlock ( "Testing Mesh generation  ..." );
  bool ok = true;
  trace.beginBlock ( "Testing Tube generation  ..." );
  std::vector<Z3i::RealPoint> aSkeleton;
  
  aSkeleton.push_back(Z3i::RealPoint(0.0, 0.0, 0.0));
  aSkeleton.push_back(Z3i::RealPoint(0.0, 0.0, 1.0));
  aSkeleton.push_back(Z3i::RealPoint(0.0, 0.0, 3.0));
  aSkeleton.push_back(Z3i::RealPoint(1.0, 2.0, 7.0));  
  aSkeleton.push_back(Z3i::RealPoint(3.0, 6.0, 10.0));  
  aSkeleton.push_back(Z3i::RealPoint(10.0, 10.0, 15.0));  
  aSkeleton.push_back(Z3i::RealPoint(12.0, 13.0, 18.0));  
  aSkeleton.push_back(Z3i::RealPoint(22.0, 15.0, 22.0));  
  aSkeleton.push_back(Z3i::RealPoint(28.0, 25.0, 32.0));  
  aSkeleton.push_back(Z3i::RealPoint(28.0, 50.0, 32.0));  
  aSkeleton.push_back(Z3i::RealPoint(28.0, 65.0, 22.0));  
  
  Mesh<Z3i::RealPoint> aMesh;
  Mesh<Z3i::RealPoint>::createTubularMesh(aMesh, aSkeleton, 2);
  for(  Mesh<Z3i::RealPoint>::VertexStorage::const_iterator it = aMesh.cVertexBegin(); 
        it !=aMesh.cVertexEnd(); it++){
    trace.info() << (*it)[0] << " " << (*it)[1] << " "<< (*it)[2] << std::endl;
  }
  std::ofstream of ("tubeMeshGeneratedFromTestMesh.off"); 
  DGtal::MeshWriter<Z3i::RealPoint>::export2OFF(of, aMesh);
  
  trace.endBlock();
  trace.info() << "Nb faces: "<< aMesh.nbFaces() << " (sould be 320)" << std::endl;
  trace.info() << "Nb vertices: "<< aMesh.nbVertex() << " (sould be 352)" << std::endl;

  bool okMeshTube1 = aMesh.nbFaces() == 320 && aMesh.nbVertex() == 352;


  ok = ok & okMeshTube1;  
  trace.endBlock();
  return ok;

}

///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

int main( int argc, char** argv )
{
  trace.beginBlock ( "Testing class Mesh" );
  trace.info() << "Args:";
  for ( int i = 0; i < argc; ++i )
    trace.info() << " " << argv[ i ];
  trace.info() << endl;

  bool res = testMesh() && testMeshGeneration(); // && ... other tests
  trace.emphase() << ( res ? "Passed." : "Error." ) << endl;
  trace.endBlock();
  return res ? 0 : 1;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
