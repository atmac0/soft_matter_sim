#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "data_types.h"

const int FIELD_WIDTH  = 100;
const int FIELD_HEIGHT = 100;

const int CELL_SIZE  = 1; //cell size in micro-meters. 


class Particle
{
private:
 
  unsigned short particle_num;      //unique number to use as a particle identifier
  Particle_type particle_type;      //type of particle, used to identify shape.
  double orientation;               //orientation measured as an angle in radians. Describes the current angle of the object, relative to the non-rotated orientation.

  unsigned int trans_velocity;      //translational velocity measured in cells per cycle
  double direction;                 //translational direction measured as an angle in radians

  double angular_velocity;          //current angular velocity, described in rads/second
  unsigned int moment_of_inertia;   //moment of inertia. dependent on mass and size of particle

  coords center_mass_coord;         //coordinate in the field of the center of mass
  coords edge_locations[];          //list of all the coordinates of the location of edges. To be used for the deletion of edged
  

public:
  void draw_square(unsigned short** field);
  void draw_edges(unsigned short** field);
  void delete_edges(unsigned short** field);
  void initialize(unsigned short par_num, Particle_type par_type, double orient, unsigned int vel, double dir, coords cm_coord);
  unsigned short get_particle_mass();
  void calc_x_y_from_angle(int &x, int &y, double angle);
};

unsigned short Particle::get_particle_mass()
{
  if(particle_type == SQUARE)
  {
    return 1;
  }
}

void Particle::initialize(unsigned short par_num, Particle_type par_type, double orient, unsigned int vel, double dir, coords cm_coord)
{
  particle_num       = par_num;
  particle_type      = par_type;
  orientation        = orient;
  trans_velocity           = vel;
  direction          = dir;
  center_mass_coord = cm_coord;
}  

/* The un-rotated orientation of the square is assumed to be as if the top and bottom of the square are parallel with the top and bottom
of the field. A line parallel to these references, pointing to the right, is used as the reference for the x-axis.

........................
.....|-----------|......
.....|...........|......
.....|.....---------->x.
.....|...........|......
.....|-----------|......
........................

The drawing begins at the bottom right corner. This location is found by taking the current orientation and subtracting pi/4 (45 degrees). Using c^2=a^2+b^2,
the x and y components of this vector are determined. This x and y component are relative to the position of the center of mass, and will be added to the coordinates
of the center of mass before being placed in the field. A value corresponding to the particle number is placed in the field. The angle value is then incremented, and the
process repeats for the whole edge, ending at the initial orientation, plus pi/4.



*/
void Particle::draw_square(unsigned short** field)
{

  double start_angle = orientation - M_PI/4;
  double end_angle   = orientation - M_PI/4;
  
  double angle = orientation - M_PI/4;
 
}

void Particle::calc_x_y_from_angle(int &x, int &y, double angle)
{
  
}


/* draw edges will identify the particle type, then call the appropriate function to draw the particle */
void Particle::draw_edges(unsigned short** field)
{
  if(particle_type == SQUARE)
  {
    draw_square(field);
  }
}

int main(){
 

  unsigned short field[FIELD_WIDTH][FIELD_HEIGHT];

  return 0;
}

void initialize_field(unsigned short** field)
{
  for(int i=0; i<FIELD_WIDTH; i++)
  {
    for(int j=0; j<FIELD_HEIGHT; i++)
    {
      field[i][j] = 0;
    }
  } 
}

// void array_to_png (char* filename, int size, array2 &field)
// {
//   char text[30];
//   int i,j;
//   double pixel;
//   pngwriter png(size,size,0,filename);
//   // Peak-finding for image scaling:
//   double max = 1e-32; 
//   for (i = 0; i < size; i++)
//   {
//     for (j = 0; j  max) max = real(field(i,j)*conj(field(i,j)));
//   }

//   // Image writeout.
//   for (i = 0; i < size; i++)
//   {
//     for (j = 0; j < size; j++)
//     {
//       pixel = real(field(i,j)*conj(field(i,j)));
//       png.plot(i,j,pixel/max,pixel/max,0.0);
//     }
//   }

//   sprintf(text, "%f7", max);
//   png.plot_text("/xtmp/dawes/share/pngwriter/fonts/FreeMonoBold.ttf",10,10,10,0.0,text,1.0,1.0,1.0);
//   png.close();
// }
