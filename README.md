# Mini-PWR-sim
A simple implementation of a simulation of a PWR (Pressurized Water Reactor). Everything is made from scratch. It has several key features often seen in full-fledged Monte Carlos Simulation engines! It samples the watt spectrum for the random assignment of energy for neutrons and more!

**FOR EDUCATIONAL PURPOSES ONLY**

# **How this works**
Firstly this doesn't emulate a PWR reactor completely, but rather is a stripped-down simple implementation!

[ It only simulates one neutron ]

So here is the full flow:

* Samples the watt spectrum to assign randomized energy values to the neutron.
* A user-specified initial position is assigned
* A direction vector is assigned (no user input, done by first choosing a random angle then using it's sin and cosine as x,y values)
* Collects needed Cross Section data from JSON database.
* Main loop mechanics:
  * Determines if the neutron is in:
      Fuel (1)
      Water moderator (2)
      Outside reactor (0)
  * Interaction Physics : 

      - In water:
    
        Calculates interaction distances
        Handles scattering events
        Updates energy and direction
      - In fuel:
        
        Handles fission events
        Simulation typically ends
  * Physical Processes
      - Scattering in Water
      - Uses realistic cross-sections
      - Models energy loss
      - Changes neutron direction
      - Tracks path length between collisions

  * Checks for:
      - Fuel rod entry
      - Moderator boundary crossing
      - System leakage

# FEATURES

- Monte Carlo Methods
    Random sampling for:
  
      - Interaction distances
  
      - Scattering angles
  
      - Fission events (based on CS data)
  
      - Physics Modeling
  
- Implements:
  
          Cross-section data
  
          Energy-dependent interactions
  
          Realistic scattering mechanics
  
- Visualization
  
    - Generates plots for:

        - Neutron trajectory
      
        - Energy degradation
  
        - System geometry
     
  IT HAS SOME ADDITIONAL FUNCTIONS THAT ARE UNUSED THAT COULD BE USED FOR FURTHER DEVELOPMENT [SOME ARE A BIT ADVANCED - FEEL FREE TO CHECK IT OUT OR CONTRIBUTE!]

  THANKS! HAVE A G'DAY! :)
