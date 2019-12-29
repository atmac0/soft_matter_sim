#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <iostream>
#include "data_types.h"
#include <string>
#include <memory>
#include <SFML/System.hpp>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>

#include "particle.h"
#include "particle_ut.h"

#include "field.h"
#include "data_types.h"

int main()
{

  Field field();
  
  particle_ut(&field);

  return 0;
}
