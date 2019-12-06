#include "particle.h"
#include "field.h"

uint16_t Particle::get_particle_mass()
{
  if(particle_type == SQUARE)
  {
    return 1;
  }
}

Particle::Particle(Field* field_, field_t par_num, Particle_type par_type, double angular_vel, double orient, double x_vel, double y_vel, coord_t cm_coord)
{
  particle_num       = par_num;
  particle_type      = par_type;

  angular_velocity   = angular_vel;
  orientation        = orient;
  
  x_velocity         = x_vel;
  y_velocity         = y_vel;

  center_mass_coord  = cm_coord;
  field              = field_;

  relative_time      = 0;
  
  if(par_type == SQUARE)
  {
    edge_locations.resize(NUM_SIDES_SQUARE*SQUARE_SIDE_LENGTH/CELL_SIZE); //resize edge location vector to hold all possible edges
    mass = MASS_SQUARE;
    moment_of_inertia = ((double)mass*pow(SQUARE_SIDE_LENGTH,2))/6; //moment of inertia is defined for a flat square as: L = (m*a^2)/6
  }
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
void Particle::add_edge_location(coord_t location)
{  
  edge_locations.push_back(location);
}

Collisions* Particle::draw_square()
{
  int32_t side_cell_count = SQUARE_SIDE_LENGTH/CELL_SIZE; /*Gives a normalized size of each cell. Aka, if the side length is 20um, and the cell size is 2um,
							the number of cells needed to draw the full edge will be 20um/2um = 10 cells */

  Collisions* collisions_p = new Collisions(SQUARE);
  coord_t coord_field;
  coord_t coord_square;
  
  for(int32_t i = -side_cell_count/2; i <= side_cell_count/2; i++)
  {
    //right side
    coord_square.x = side_cell_count/2;
    coord_square.y = i;
    coord_field = translate_to_field(coord_square);
    field->place_edge_in_field(coord_field, particle_num, collisions_p);
    add_edge_location(coord_field);

    //top side
    coord_square.x = i;
    coord_square.y = -side_cell_count/2;
    coord_field = translate_to_field(coord_square);
    field->place_edge_in_field(coord_field, particle_num, collisions_p);
    add_edge_location(coord_field);
    
    //left side
    coord_square.x = -side_cell_count/2;
    coord_square.y = i;
    coord_field = translate_to_field(coord_square);
    field->place_edge_in_field(coord_field, particle_num, collisions_p);
    add_edge_location(coord_field);
    
    //bottom side
    coord_square.x = i;
    coord_square.y = side_cell_count/2;
    coord_field = translate_to_field(coord_square);
    field->place_edge_in_field(coord_field, particle_num, collisions_p);
    add_edge_location(coord_field);
  }
  return collisions_p;
}

//Takes a coordinate of an edge point relative to the center of mass of the particle, and translates it into coordinates relative to the field
coord_t Particle::translate_to_field(coord_t edge_point)
{

  /*From a 2d matrix transform:
    x' = xcos(theta) - ysin(theta)
    y' = xsin(theta) + ycos(theta)
    where x and y are the initial coordinates, and x' and y' are the rotated coordinates */
  int32_t x_prime = edge_point.x*cos(orientation) - edge_point.y*sin(orientation); 
  int32_t y_prime = edge_point.x*sin(orientation) + edge_point.y*cos(orientation);

  //Translate from coordinates relative to the center of mass, to relative to the field
  edge_point.x = center_mass_coord.x + x_prime;
  edge_point.y = center_mass_coord.y + y_prime;	 

  //Check if edge is inside of the bounds of the field. If not, have it wrap around.
  if(edge_point.x < 0)
  {
    edge_point.x = FIELD_WIDTH + edge_point.x;
  }
  else
  {
    edge_point.x = edge_point.x%FIELD_WIDTH;
  }

  if(edge_point.y < 0)
  {
    edge_point.y = FIELD_HEIGHT + edge_point.y;
  }
  else
  {
    edge_point.y = edge_point.y%FIELD_HEIGHT;
  }
  
  
  return edge_point;
}

/* draw edges will identify the particle type, then call the appropriate function to draw the particle */
Collisions* Particle::draw_edges()
{
  if(particle_type == SQUARE)
  {
    return draw_square();
  }
}

void Particle::delete_edges()
{
  for(std::vector<coord_t>::iterator it = edge_locations.begin(); it != edge_locations.end(); it++)
  {
    if(field->get_particle_at(*it) == particle_num)
    {
      field->clear_particle_at(*it);
    }
  }

  edge_locations.clear();  
}

void Particle::propagate()
{

  double time = 0; //time spent propogating

  double x_cell_velocity = abs(x_velocity)/CELL_SIZE; //convert x velocity to cells per micro-second
  double y_cell_velocity = abs(y_velocity)/CELL_SIZE; //convert y velocity to cells per micro-second

  double x_cell_time     = 1/x_cell_velocity; //time spent per cell movement
  double y_cell_time     = 1/y_cell_velocity;

  double x_distance_to_travel = y_cell_time * abs(x_velocity);
  double y_distance_to_travel;

  uint32_t x_cells_to_travel;
  uint32_t y_cells_to_travel;
  
  double x_cell_int_as_double; //used to store the integer portion of the x granular position
  double y_cell_int_as_double; //used to store the integer portion of the y granular position

  double x_granular_decimal; //used to store the portion of the x granular position after the decimal
  double y_granular_decimal; //used to store the portion of the y granular position after the decimal
  
  Collisions* collisions;

  /*propogate the particle for the TIME_INCREMENT amout of time. Curret time is kept track of relative to the change in the y coordinate.
   For every cell moved in the y direction, move x_y_ratio in the x direction. The movement in the y direction will occur first, followed
   by a redraw. For every whole integer in x_y_ratio, the particle will be propagated, with a redraw at each cell, in the x direction that
  number of cells. The rollover will be stored in the granular position variable of the particle.*/
  while(time < TIME_INCREMENT)
  {
    /*Check if the particle will travel atleast 1 cell in the y direction in the time remaining. If it does, proceed by incrementing
      the y direction by one, and the x direction by the number it should travel, determined by the x_y_ratio.*/
    
    if((time + y_cell_time) <= TIME_INCREMENT)
    {

      rotate_particle(y_cell_time);
      
      collisions = translate_y_by_1();
      resolve_collisions(collisions);

      /*being a little tricky here. translate_x_by_granular will return either 0 or 1 if the cell overflows. the modf will take the integer
	portion of the double x_distance_to_travel, which will be assigned to x_cells_to_travel. If the decimal remainder overflows the the cell,
	the number of cells to travel will then have 1 added to it.*/
      
      x_cells_to_travel = 0;
      x_cells_to_travel += translate_x_by_granular(modf(x_distance_to_travel, &x_cell_int_as_double));
      x_cells_to_travel += (uint32_t)x_cell_int_as_double;

      for(uint32_t i=0; i<x_cells_to_travel; i++)
      {
      	collisions = translate_x_by_1();
      	resolve_collisions(collisions);
      }
      
      time += y_cell_time;
    }
    /*If incrementing once in the y direction would go over the time increment limit, add the granular amount. If the granular amount goes over one,
      propogate in the y direction. Then propogate x by the correct proportion.*/
    else if((time + y_cell_time) > TIME_INCREMENT)
    {
      double time_remaining = ((double)TIME_INCREMENT) - time;

      rotate_particle(time_remaining);
      
      y_distance_to_travel = abs(y_velocity) * time_remaining;
      x_distance_to_travel = abs(x_velocity) * time_remaining;

      if(translate_y_by_granular(y_distance_to_travel) == 1)
      {
    	collisions = translate_y_by_1();
	resolve_collisions(collisions);
      }
      else
      {
	//This is only called if the y velocity is very small. Without this, the particle will not be drawn if it is stationary.
	delete_edges();
	resolve_collisions(draw_edges());
      }

      x_cells_to_travel = 0;
      x_cells_to_travel += translate_x_by_granular(modf(x_distance_to_travel, &x_cell_int_as_double));
      x_cells_to_travel += (uint32_t)x_cell_int_as_double;

      for(uint32_t i = 0; i < x_cells_to_travel; i++)
      {
    	collisions = translate_x_by_1();
    	resolve_collisions(collisions);
      }
      time = TIME_INCREMENT;
    }
  }
}

//add the rotational displacement (angular velocity * time span) to the orientation. If the orientation goes over 2*PI, take the remainder.
void Particle::rotate_particle(double time_span)
{
  orientation += angular_velocity * time_span;
  orientation = fmod(orientation, M_PI*2);
}

/*Translate the particle by 1 in the direction appropriate to the current
velocity of the particle in the x direction. This fuction deletes the current edges, propagates
the particle in the appropriate direction by 1, then redraws the edges. This
function returns the collisions that occured from propogation.

If the particle goes outside of the field, it will wrap around, meaning that if the particle exits
the bottom, it will reenter thetop, or if it exits the left, it will reenter the right.*/
Collisions* Particle::translate_x_by_1()
{

  delete_edges();
  
  if(x_velocity > 0)
  {
    center_mass_coord.x++;
    center_mass_coord.x = center_mass_coord.x%FIELD_WIDTH;
  }
  else
  {
    center_mass_coord.x--;
    center_mass_coord.x = center_mass_coord.x%FIELD_WIDTH;
  }

  return draw_edges();
}

/*Translate the particle by 1 in the direction appropriate to the current
velocity of the particle in the y direction. This fuction deletes the current edges, propagates
the particle in the appropriate direction by 1, then redraws the edges. This
function returns the collisions that occured from propogation.

If the particle goes outside of the field, it will wrap around, meaning that if the particle exits
the bottom, it will reenter thetop, or if it exits the left, it will reenter the right.*/
Collisions* Particle::translate_y_by_1()
{
  
  delete_edges();
  
  if(y_velocity > 0)
  {
    center_mass_coord.y++;
    center_mass_coord.y = center_mass_coord.y%FIELD_HEIGHT;
  }
  else
  {
    center_mass_coord.y--;
    center_mass_coord.y = center_mass_coord.y%FIELD_HEIGHT;
  }

  return draw_edges();
}

/*Translates a particle by a granular distance. If the granular postion
  does not exceed the bounds of the cell (is not greater than 1 or less than 0), 
the function will return 0, meaning 0 cells should be moved. If it does exceed those 
bounds, the fuction will return 1, meaning the center of mass postion needs to be moved
by 1.*/
uint32_t Particle::translate_y_by_granular(double granularity)
{
  if(y_velocity > 0)
  {
    y_position_granular += granularity;
  }
  else
  {
    y_position_granular -= granularity;
  }

  if(y_position_granular < 0)
  {
    y_position_granular = 1 + y_position_granular;
    return 1;
  }
  if(y_position_granular > 1)
  {
    y_position_granular = y_position_granular - 1;
    return 1;
  }

  return 0;
}

uint32_t Particle::translate_x_by_granular(double granularity)
{
  if(x_velocity > 0)
  {
    x_position_granular += granularity;
  }
  else
  {
    x_position_granular -= granularity;
  }

  if(x_position_granular < 0)
  {
    x_position_granular = 1 + x_position_granular;
    return 1;
  }
  if(x_position_granular > 1)
  {
    x_position_granular = x_position_granular - 1;
    return 1;
  }

  return 0;
}

void Particle::resolve_collisions(Collisions* collisions)
{
    //dp_t: change in translational momentum; dp_l: change in angular momentum. 1 is the current particle in iteration.
    // 2 is the particle in which the collision happened.
  momentum_t delta_translational[NUM_PARTICLES] = {0};
  double delta_angular[NUM_PARTICLES] = {0};
  momentum_t p1, p2;

  coord_t collision_location;
  field_t particle_1_num, particle_2_num;
  coord_t particle_1_cm_coord, particle_2_cm_coord;

  //Iterate through all the collisions that occured during the most recent drawing
  for(uint32_t i = 0; i < collisions->counts(); i++)
  {   
    std::cout << "Collision at (" << collisions->get_collision_at(i).location.x << "," << collisions->get_collision_at(i).location.y << ") between particle " << collisions->get_collision_at(i).particle1 << " and particle " << collisions->get_collision_at(i).particle2 << "!\n";

    //Get the location of the collision and the identifying number of the particles that interacted
    collision_location = collisions->get_collision_at(i).location;
    particle_1_num             = collisions->get_collision_at(i).particle1;
    particle_2_num             = collisions->get_collision_at(i).particle2;

    //Get the center of mass positions of each particle
    particle_1_cm_coord        = particles[particle_1_num]->center_mass_coord;
    particle_2_cm_coord        = particles[particle_2_num]->center_mass_coord;

    //Find the linear momentum of the edge of each particle at the collision location
    p1 = find_linear_momentum_at(collision_location, particle_1_num);
    p2 = find_linear_momentum_at(collision_location, particle_2_num);

    //Find the transfter of momentum due to the individual edge points that collided. The final change will be the sum of the change due to
    //every colliding edge point.
    std::cout << "Change in momentum of particle " << particle_2_num << " due to collision at ("
	      << collision_location.x << "," << collision_location.y << "):\n";
    find_change_in_momentum(p1, particle_2_cm_coord, collision_location, &delta_translational[particle_2_num], &delta_angular[particle_2_num]);
    std::cout << "Change in momentum of particle " << particle_1_num << " due to collision at ("
	      << collision_location.x << "," << collision_location.y << "):\n";
    find_change_in_momentum(p2, particle_1_cm_coord, collision_location, &delta_translational[particle_1_num], &delta_angular[particle_1_num]);
  }

  //This is ~dirty~. I had to use a list of length NUM_PARTICLES because I needed some way to keep track of which particles were interacted with.
  //Most of the elements of the list are going to be empty, however I still iterate through all of them. Implementing some way to keep track of
  //just the particles that were interacted with would be good.
  for(uint32_t i = 0; i < NUM_PARTICLES; i++)
  {
    Particle* particle = particles[i];
    //To find the velocity, p=mv, where v=p/m. This is why you must divide by momentum by the mass.
    if(collisions->counts() != 0)
    {
      double average_dx_trans_vel = (delta_translational[i].x/particle->get_mass())/collisions->counts();
      double average_dy_trans_vel = (delta_translational[i].y/particle->get_mass())/collisions->counts();
      particle->increment_x_velocity(average_dx_trans_vel);
      particle->increment_y_velocity(average_dy_trans_vel);
    }
  }
}

/*The change in momentum has two parts, the change in translational momentum, and the change in angular momentum. To find the translational change, 
find the projection of the linear momentum vector with the vector pointing from the colliding edge to the center of mass. The x and y components of this
vector will be the */
void Particle::find_change_in_momentum(momentum_t linear_momentum, coord_t cm_coord, coord_t coll_location, momentum_t* dp_t, double* dp_l)
{
  uint32_t x_diff    = cm_coord.x - coll_location.x;
  uint32_t y_diff    = cm_coord.y - coll_location.y;
  double   r_mag = sqrt(pow(x_diff,2) + pow(y_diff,2));
  double lin_mom_mag = sqrt(pow(linear_momentum.x,2) + pow(linear_momentum.y,2));  
  
  dp_t->x += (linear_momentum.x * x_diff)/r_mag;
  dp_t->y += (linear_momentum.y * y_diff)/r_mag;
 
  double dp_lx = (linear_momentum.x *  y_diff)/r_mag;
  double dp_ly = (linear_momentum.y * -x_diff)/r_mag;
  
  *dp_l += sqrt( pow(dp_lx,2) + pow(dp_ly,2) );

  std::cout << "Change in linear momentum: dx = " << (linear_momentum.x * x_diff)/r_mag << ", dy = " << (linear_momentum.y * y_diff)/r_mag << "\n";
  std::cout << "Change in angular momntum: dl = " << sqrt( pow(dp_lx,2) + pow(dp_ly,2) ) << "\n";
}

/*Find the linear momentum at a point on the edge of a particle.
Notation:
L: angular momentum
v: translational velocity
r: disance from point of rotation (center of mass)
w: angular velocity
m: mass of particle
I: moment of inertia

L = rp; L=Iw; => p = Iw/r

These relations are used to find the magnitude of the momentum at a point on the edge. The components
are then found using
*/
momentum_t Particle::find_linear_momentum_at(coord_t point, field_t particle_num)
{

  Particle* particle = particles[particle_num];
  momentum_t linear_momentum;

  //Find the components due to translational momentum
  linear_momentum.x = particle->get_x_velocity()*particle->get_mass();
  linear_momentum.y = particle->get_y_velocity()*particle->get_mass();

  //Find the components due to angular momentum
  /*  To find the linear momentum due to the angular velocity, imagine you draw a circle about the center of mass where
  the edge intersects with the point at which a collision happened. Draw a line from the center of mass to this point.
  The direction of the momentum is the negative inverse of the slope of this line. However, the direction is opposite if 
  the particle is rotating counter clockwise.

  ................+.......
  .....|-----------*+..... The stars represent the r vector. The + represent the perpendicular momentum vector.
  .....|........*..|..+....
  .....|.....*--------->x.
  .....|...........|......
  .....|-----------|......
  ........................
*/

  uint32_t x_diff = point.x - particle->get_center_mass_coord().x;
  uint32_t y_diff = point.y - particle->get_center_mass_coord().y;

  double w        = particle->get_angular_velocity();
  double I        = particle->get_moment_of_inertia();
  double r        = sqrt(pow(x_diff,2) + pow(y_diff,2));
  
  double slope;    
  double p_slope;   
  
  double p_magnitude = I*w/r;

  //Check if the difference in the x or y is zero. If it is, all of the angular momentum is linear in the corresponding x or y direction.
  if(x_diff == 0)
  {
    linear_momentum.x += p_magnitude;
    return linear_momentum;
  }
  if(y_diff == 0)
  {
    linear_momentum.y += p_magnitude;
    return linear_momentum;
  }

  slope = y_diff/x_diff;
  p_slope = 1/slope;

  double theta = atan(p_slope);
  linear_momentum.x += p_magnitude*sin(theta);
  linear_momentum.y += p_magnitude*cos(theta);
  return linear_momentum;
}

uint32_t Particle::get_mass()
{
  return mass;
}

void Particle::increment_y_velocity(double dy)
{
  y_velocity += dy;
}

void Particle::increment_x_velocity(double dx)
{
  x_velocity += dx;
}

void Particle::increment_angular_velocity(double dw)
{
  angular_velocity += dw;
}

double Particle::get_y_velocity()
{
  return y_velocity;
}

double Particle::get_x_velocity()
{
  return x_velocity;
}

coord_t Particle::get_center_mass_coord()
{
  return center_mass_coord;
}

double Particle::get_angular_velocity()
{
  return angular_velocity;
}

void Particle::set_particles_array(Particle** par_arr)
{
  particles = par_arr;
}

double Particle::get_moment_of_inertia()
{
  return moment_of_inertia;
}
