#include "particle.h"
#include "field.h"

uint16_t Particle::get_particle_mass()
{
  if(particle_type == SQUARE)
  {
    return 1;
  }
}

Particle::Particle(Field* field_, uint16_t par_num, Particle_type par_type, double orient, uint32_t trans_vel, double angular_vel, double dir, coord_t cm_coord)
{
  particle_num       = par_num;
  particle_type      = par_type;
  orientation        = orient;
  trans_velocity     = trans_vel;
  angular_velocity   = angular_vel;
  direction          = dir;
  center_mass_coord  = cm_coord;
  field              = field_;
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

The drawing begins at the bottom right corner. This location is found by converting from polar coordinates (using the angle of the orientation),
to cartesian coordinates (using the field as a reference frame), where:

x = rcos(theta)
y = rsin(theta)

where r is calulated by incremented values of an x and y value that are stored relative to the square, using the pythagorean theorem.

*/
Collisions Particle::draw_square()
{
  
  int32_t side_cell_count = SQUARE_SIDE_LENGTH/CELL_SIZE; /*Gives a normalized size of each cell. Aka, if the side length is 20um, and the cell size is 2um,
							the number of cells needed to draw the full edge will be 20um/2um = 10 cells */


  Collisions collisions(SQUARE);
  coord_t coord_field;
  coord_t coord_square;
  
  for(int32_t i=-side_cell_count/2; i<=side_cell_count/2; i++)
  {
    //right side
    coord_square.x = side_cell_count/2;
    coord_square.y = i;
    coord_field = translate_to_field(coord_square);
    field->place_edge_in_field(coord_field, particle_num, &collisions);

    //top side
    coord_square.x = i;
    coord_square.y = -side_cell_count/2;
    coord_field = translate_to_field(coord_square);
    field->place_edge_in_field(coord_field, particle_num, &collisions);
    
    //left side
    coord_square.x = -side_cell_count/2;
    coord_square.y = i;
    coord_field = translate_to_field(coord_square);
    field->place_edge_in_field(coord_field, particle_num, &collisions);
    
    //bottom side
    coord_square.x = i;
    coord_square.y = side_cell_count/2;
    coord_field = translate_to_field(coord_square);
    field->place_edge_in_field(coord_field, particle_num, &collisions);
  }
}

//Takes a coordinate of an edge point relative to the center of mass of the particle, and translates it into coordinates relative to the field
coord_t Particle::translate_to_field(coord_t edge_point)
{

  /*From a 2d matrix transform:
    x' = xcos(theta) - ysin(theta)
    y' = xsin(theta) + ycos(theta)
    where x and y are the initial coordinates, and x' and y' are the rotated coordinates */
  int32_t x_prime = edge_point.x*cos(-orientation) - edge_point.y*sin(-orientation); 
  int32_t y_prime = edge_point.x*sin(-orientation) + edge_point.y*cos(-orientation);

  //Translate from coordinates relative to the center of mass, to relative to the field
  edge_point.x = center_mass_coord.x + x_prime;
  edge_point.y = center_mass_coord.y + y_prime;  
  
  return edge_point;
}

/* draw edges will identify the particle type, then call the appropriate function to draw the particle */
Collisions Particle::draw_edges()
{
  if(particle_type == SQUARE)
  {
    return draw_square();
  }
}
