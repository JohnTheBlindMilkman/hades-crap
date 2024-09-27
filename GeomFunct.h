#include "TMath.h"
#include "hgeomvector.h"
#include "hphysicsconstants.h"


#include <iostream>
#include <stdio.h>
#include <math.h>

using namespace std;


void CalcSegVector(Double_t z, Double_t rho, Double_t phi, Double_t theta, HGeomVector &base, HGeomVector &dir)
{
    static const Double_t pi2=TMath::PiOver2();
    base.setXYZ(rho*cos(phi+pi2), rho*sin(phi+pi2), z);
    dir.setXYZ(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
};


Double_t calcRMS(const Double_t* valArr, Double_t Mean,Int_t valNum)
{
    // Calculates RMS of valNum numbers in valArr using the Mean of the values
    Double_t RMS = 0.0;
    Double_t sum = 0.0;

    for(Int_t i = 0; i < valNum; i++)
    {
        sum += pow(valArr[i]-Mean,2.0);
    }
    RMS = sqrt(sum/valNum);
    return RMS;
}

//versions from functions_TS_Alex.h:
Double_t calcDeterminant(HGeomVector& v1, HGeomVector& v2, HGeomVector& v3)
{
  // calculating the Determinant of a 3 x 3 Matrix
  // with the column vectors [v1, v2, v3]
  // using the RULE of SARRUS
  //
  // | v1(0)   v2(0)   v3(0) |      | v1(0) v2(0) v3(0)|v1(0) v2(0)  .
  // |                       |      |  \\     \\     X   |  /     /    .
  // |                       |      |   \\     \\   / \\  | /     /     .
  // |                       |      |    \\     \\ /   \\ |/     /      .
  // |                       |      |     \\     X     \\/     /       .
  // |                       |      |      \\   / \\    /\\    /        . 
  // |                       |      |       \\ /   \\  / |\\  /         .
  // | v1(1)   v2(1)   v3(1) |   =  | v1(1) v2(1) v3(1)|v1(1) v2(1)  .
  // |                       |      |       / \\    /\\  | /\\          .
  // |                       |      |      /   \\  /  \\ |/  \\         .
  // |                       |      |     /     \\/    \\/    \\        .
  // |                       |      |    /      /\\    /\\     \\       .
  // |                       |      |   /      /  \\  / |\\     \\      . 
  // |                       |      |  /      /    \\/  | \\     \\     .
  // | v1(2)   v2(2)   v3(2) |      | v1(2) v2(2) v3(2)| v1(2) v2(2) . 
  //                                 /      /     /  \\     \\     \\   .
  //                                                               
  //                                -      -     -    +     +     +  .

  return ( v1(0) * v2(1) * v3(2)
           + v2(0) * v3(1) * v1(2)
           + v3(0) * v1(1) * v2(2)
           - v3(0) * v2(1) * v1(2)
           - v1(0) * v3(1) * v2(2)
           - v2(0) * v1(1) * v3(2));
}


/*
Double_t calculateMinimumDistanceStraightToPoint(Double_t &baseX, Double_t &baseY, Double_t &baseZ, Double_t &dirX, Double_t &dirY, Double_t &dirZ, Double_t &pointX, Double_t &pointY, Double_t &pointZ)

{
    // calculates the minimum distance of a point to a straight given as parametric straight x = base + n * dir
    HGeomVector base, dir, point;
    base.setXYZ(baseX,baseY,baseZ);
    dir.setXYZ(dirX,dirY,dirZ);
    point.setXYZ(pointX,pointY,pointZ);

    if (!(dir.length()>0))
    {
        return -1000000.;
    }
 
    HGeomVector diff = base-point;

    HGeomVector cross = dir.vectorProduct(diff);
 
    return cross.length()/dir.length();
}
*/

//  other version:

Double_t calculateMinimumDistanceStraightToPoint(HGeomVector &base, HGeomVector &dir,HGeomVector &point)

{
   //calculates the minimum distance of a point to a straight given as parametric straight x = base + n * dir

  if (!(dir.length()>0))
    {
      return -1000000.;
    }
 
  HGeomVector diff = base-point;

  HGeomVector cross = dir.vectorProduct(diff);
 
  return cross.length()/dir.length();
}



//--------------------------------------------------------------------------------------------------------
Double_t calculateMinimumDistance(HGeomVector &base1, HGeomVector &dir1, HGeomVector &base2, HGeomVector &dir2)
{
  // calculates the minimum distance of two tracks given as parametric straights x = base + n * dir

  HGeomVector cross = dir1.vectorProduct(dir2);

  HGeomVector ab = base1 - base2;

  if ( !( fabs(cross.length())>0.)) // dir1 || dir2
  {
      return calculateMinimumDistanceStraightToPoint(base1, dir1, base2);

      //Double_t base1x = base1.getX();
      //Double_t base1y = base1.getY();
      //Double_t base1z = base1.getZ();
      //Double_t dir1x  = dir1.getX();
      //Double_t dir1y  = dir1.getY();
      //Double_t dir1z  = dir1.getZ();
      //Double_t base2x = base2.getX();
      //Double_t base2y = base2.getY();
      //Double_t base2z = base2.getZ();

     // return calculateMinimumDistanceStraightToPoint(base1x,base1y,base1z,dir1x,dir1y,dir1z,base2x,base2y,base2z);
     // return calculateMinimumDistanceStraightToPoint(base1.getX(),base1.getY(),base1.getZ(),dir1.getX(),dir1.getY(),dir1.getZ(),base2.getX(),base2.getY(),base2.getZ());

      //cout<<"inside IF 2"<<endl;
  }
  return fabs(ab.scalarProduct(cross)/cross.length());
}


HGeomVector calculatePointOfClosestApproach(HGeomVector &base1, HGeomVector &dir1,
                                                                    HGeomVector &base2, HGeomVector &dir2)
{
  //  calculating point of closest approach
  //       
  //        from the equations of the straight lines of g and h
  //        g: x1 = base1 + l * dir1
  //        h: x2 = base2 + m * dir2
  //       
  //        you can construct the following planes:
  //       
  //        E1: e1 = base1  +  a * dir1  +  b * (dir1 x dir2)
  //        E2: e2 = base2  +  s * dir2  +  t * (dir1 x dir2)
  //       
  //        now the intersection point of E1 with g2 = {P1}
  //        and the intersection point of E2 with g1 = {P2}
  //       
  //        form the base points of the perpendicular to both straight lines.
  //       
  //        The point of closest approach is the middle point between P1 and P2:
  //       
  //        vertex = (p2 - p1)/2
  //
  //        E1 ^ g2:
  //
  //           e1 = x2
  //    -->    base1  +  a * dir1  +  b * (dir1 x dir2) = base2 + m * dir2
  //    -->    base1 - base2 = m * dir2  -  a * dir1  -  b * (dir1 x dir2)      
  //                                          (m)
  //    -->    [ dir2, -dir1, -(dir1 x dir2)] (a) = base1 - base2       
  //                                          (b)
  //          
  //           using CRAMER's RULE you can find the solution for m (a,b, not used)
  //          
  //           using the rules for converting determinants:
  //          
  //           D12 = det [dir2, -dir1, -(dir1 x dir2)]
  //               = det [dir2,  dir1,  (dir1 x dir2)]
  //          
  //           Dm  = det [base1 - base2, -dir1, -(dir1 x dir2)]
  //               = det [base1 - base2,  dir1,  (dir1 x dir2)]
  // 
  //            m  = Dm/D12
  //          
  //           P1: p1 = x2(m)
  //                  = base2 + Dm/D12 * dir2
  //
  //        E2 ^ g1:
  //
  //           e2 = x1
  //    -->    base2  +  s * dir2  +  t * (dir1 x dir2) = base1 + l * dir1
  //    -->    base2 - base1 = l * dir1  -  s * dir2  -  t * (dir1 x dir2)      
  //                                          (l)
  //    -->    [ dir1, -dir2, -(dir1 x dir2)] (s) = base2 - base1       
  //                                          (t)
  //          
  //           again using CRAMER's RULE you can find the solution for m (a,b, not used)
  //          
  //           using the rules for converting determinants:
  //          
  //           D21 =  det [dir1, -dir2, -(dir1 x dir2)]
  //               =  det [dir1,  dir2,  (dir1 x dir2)]
  //               = -det [dir2,  dir1,  (dir1 x dir2)]
  //               = -D12
  //          
  //           Dl  =  det [base2 - base1, -dir2, -(dir1 x dir2)]
  //               =  det [base2 - base1,  dir1,  (dir1 x dir2)]
  //               = -det [base1 - base2,  dir1,  (dir1 x dir2)]
  //
  //            l  =   Dl/D21
  //               = - Dl/D12
  //          
  //           P2: p2 = x1(m)
  //                  = base1 - Dl/D12 * dir1
  //          
  //          
  //           vertex = p1 + 1/2 * (p2 - p1)
  //                  = 1/2 * (p2 + p1)
  //                  = 1/2 *( (base1 + base2) +  1/D12 * ( Dm * dir2 - Dl * dir1) )
  //                     

  HGeomVector cross = dir1.vectorProduct(dir2); // cross product: dir1 x dir2

  // straight lines are either skew or have a cross point
             
  HGeomVector diff = base1;
  diff-=base2; // Difference of two base vectors base1 - base2
               
  Double_t D;
  D =  calcDeterminant(dir2, dir1 ,cross);

  if (!(fabs(D) > 0.))
    {
      ::Warning(":calculatePointOfClosestApproach","Dirs and cross-product are lin. dependent: returning default Vertex (-20000,-20000,-20000)");

      return HGeomVector(-20000.,-20000.,-20000.);
    }

  Double_t Dm =  calcDeterminant(diff , dir1, cross);
  Double_t Dl = -calcDeterminant(diff , dir2, cross);

  HGeomVector vertex;
  HGeomVector dm;
  HGeomVector dl;

  dm = dir2;
  dm *= Dm;

  dl = dir1;
  dl *= Dl;

  vertex = dm - dl;

  vertex *= ((1.)/D);

  vertex+=base1;
  vertex+=base2;
  vertex*=0.5;

  return HGeomVector(vertex);
}

//--------------------------------------------------------------------------------------------------------

HGeomVector calculateCrossPoint(HGeomVector &base1, HGeomVector &dir1,
                                                        HGeomVector &base2, HGeomVector &dir2)
{
  // calculating cross point
  // taking all three equations into account solving the overdetermined set of lin. equations
  // of
  // base1 + l * dir2 =  base1 + m * dir2
  //
  // set of lin. equations:
  // 
  //   base1(0) + l * dir1(0) = base2(0) + m * dir2(0)
  //   base1(1) + l * dir1(1) = base2(1) + m * dir2(1)
  //   base1(2) + l * dir1(2) = base2(2) + m * dir2(2) this line is ignored
  //
  //   written in matrix form
  //
  //        l \\\
  //   M * |   | = base2 - base1
  //       \\ m /
  //
  //   M is a 3x2 matrix
  //    
  // to solve multiply the equation by the transposed Matrix of M from the left: M
  //    
  //  T      /  l \\                                                               .
  // M * M * |    | = M  * (base2 - base1)
  //         \\ -m /
  // MIND THE '-' of m
  //
  //     / dir1(0) dir2(0) \\                                                      .
  //     |                 |    T   / dir1(0) dir1(1) dir1(2) \\                   .
  // M = | dir1(1) dir2(1) |,  M  = |                         |
  //     |                 |        \\ dir2(0) dir2(1) dir2(2) /                   .
  //     \\ dir1(2) dir2(2) /                                   
  //
  //  T      / (dir1(0)*dir1(0) + dir1(1)*dir1(1) + dir1(2)*dir1(2))   (dir1(0)*dir2(0) + dir1(1)*dir2(1) + dir1(2)*dir2(2))  \\ .

  // M * M = |                                                                                                                |

  //         \\ (dir1(0)*dir2(0) + dir1(1)*dir2(1) + dir1(2)*dir2(2))   (dir2(0)*dir2(0) + dir2(1)*dir2(1) + dir2(2)*dir2(2))  /                       

  //
  //  T       / d1d1 d1d2 \\                           .
  // M  * M = |           |
  //          \\ d1d2 d2d2 /
  //
  // diff = base2 - base1
  //
  //  T           /  (dir1(0)*diff(0) + dir1(1)*diff(1) + dir1(2)*diff(2)) \\         .
  // M  * diff =  |                                                        |
  //              \\  (dir2(0)*diff(0) + dir2(1)*diff(1) + dir2(2)*diff(2)) /
  //
  //  T           /  d1diff  \\                                          .
  // M  * diff =  |          |
  //              \\  d2diff  /
  //
  // now the new Matrix set is to be solved by CRAMER'S Rule:
  //
  // / d1d1 d1d2 \\   /  l \\   /  d1diff \\                   .
  // |           | * |    | = |          |
  // \\ d1d2 d2d2 /   \\ -m /   \\  d2diff /
  //
  //     | d1d1 d1d2 |
  // D = |           | = d1d1*d2d2 - d1d2*d1d2;
  //     | d1d2 d2d2 |
  //
  //     | d1diff d1d2 |
  // Dl= |              | = d1diff*d2d2 - d1d2*d2diff;
  //     | d2diff d2d2 |             
  //
  // l = Dl/D = l_cross
  //
  // vertex = base1 + l_cross * dir1
  //

  Double_t d1d1 = dir1(0)*dir1(0) + dir1(1)*dir1(1) + dir1(2)*dir1(2);
  Double_t d2d2 = dir2(0)*dir2(0) + dir2(1)*dir2(1) + dir2(2)*dir2(2);
  Double_t d1d2 = dir1(0)*dir2(0) + dir1(1)*dir2(1) + dir1(2)*dir2(2);
 
  Double_t D = d1d1*d2d2 - (d1d2*d1d2);
 
  if (!(fabs(D) > 0.))
    {
      ::Warning("calculateCrossPoint","Error while calculating cross point ... eqns are lin. dependent:returning default Vertex (-20000,-20000,-20000)");

      return HGeomVector(-20000.,-20000.,-20000.);
    }

  Double_t d1diff = dir1(0)*(base2(0)-base1(0))+dir1(1)*(base2(1)-base1(1))+dir1(2)*(base2(2)-base1(2));
  Double_t d2diff = dir2(0)*(base2(0)-base1(0))+dir2(1)*(base2(1)-base1(1))+dir2(2)*(base2(2)-base1(2));

  Double_t Dlambda = d1diff*d2d2-d1d2*d2diff;
 
  Double_t lambda = Dlambda/D;
 
  HGeomVector vertex;
  vertex += dir1;
  vertex *= lambda;
  vertex += base1;

  cout << "Cross point calculated" << endl;
  return HGeomVector(vertex);

 // return HGeomVector(-20000.,-20000.,-20000.);
}

HGeomVector calcVertexAnalytical(HGeomVector &base1, HGeomVector &dir1,
                                                         HGeomVector &base2, HGeomVector &dir2)
{
  // Calculates the Vertex of two straight lines defined by the vectors base and dir
  //
  //      g: x1 = base1 + l * dir1
  //      h: x2 = base2 + m * dir2 , where l,m are real numbers
  //                                   h
  // 1. are g and h
  //       parallel / identical, i.e. are dir1 and dir2 linear dependent?
  //      
  //                                        /-                              
  //                                        |
  //                                        |   = 0    linear dependent, no unique solution, returning dummy 
  //      => cross product : dir1 x dir2 = -| 
  //                                        |  != 0    linear independent
  //                                        |
  //                                        \\-        
  //
  // 2. are g and h
  //       skew or do they have a crossing point, i.e are dir1, dir2 and (base1 - base2) linear dependent ?
  //
  //                                                    /-                              
  //                                                    |
  //                                                    |   = 0    linear dependent
  //                                                    |          g and h are intersecting
  //                                                    |          calculating vertex as point of intersection
  //                                                    |
  //    => determinant: det[ dir1, dir2, base1-base2]= -|
  //                                                    |  != 0    linear independent
  //                                                    |          g and h are skew
  //                                                    |          calulating vertex as point of closest approach
  //                                                    |
  //                                                    \\-        
  // 
  // 3.
  //    (a) calculating intersection point
  //    (b) calculating point of closest approach


  // 1. exists a unique solution ?

  if ((dir1.vectorProduct(dir2)).length()> 0.) // dir1 and dir2 linear independent
    {
      // straight lines are either skew or have a cross point

      HGeomVector diff = base1;
      diff-=base2; // Difference of two base vectors base1 - base2
     
      // 2. skew or intersecting ?
       
      if (fabs(calcDeterminant(dir2, dir1 ,diff))>0.)
        {
          // 3. (b) skew
          return HGeomVector(calculatePointOfClosestApproach(base1, dir1, base2, dir2));
        }
      else
        {
          // 3. (a) intersection
          return HGeomVector(calculateCrossPoint(base1 ,dir1, base2 ,dir2));
        }
    }
  else
    {
      // dir1 and dir2 linear dependent -> g1 and g2 identical or parallel
      return HGeomVector(-10000000.,-10000000.,-10000000.);
    }
  return HGeomVector(-10000000.,-10000000.,-10000000.);
} 