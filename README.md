# ReactionDiffusionSimulator
Simulates Reaction Diffusion models
This tool simulates a number of reaction-diffusion systems and produces Turing patterns. Optionally, the images may be saved sequentially and can output images to string together in an animation.

##Examples 

#### Example usage:

    python ReactionDiffusion.py -o my_simulation --moviemode -n 10000

This will output files with the prefix "OutputPrefix" into a directory of the same name, and as
the moviemode flag is set, it will store 250 images for animation. User will be prompted for model type

#### Example usage 2:

    python ReactionDiffusion.py -o my_simulation2 -m GM -n 5000

Instead uses the Gierer-Meinhardt activator-inhibitor model (-gm).

## Supported models
Currently supports the following reaction-diffusion systems

### Gray-Scott
####Solitons
<img src="https://github.com/jluebeck/ReactionDiffusionSimulator/blob/master/images/solitons.png" width="640">

####Coral

####Maze

####Waves

####Flicker

####Worms

### FitzHugh-Nagumo

### Gierer-Meinhardt

