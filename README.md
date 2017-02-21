# Why I did this fork:
For my package [climex](https://github.com/theGreatWhiteShark/climex) I needed an optimization function, that provides me with ALL the updates it is performing in order to do an animation.

Since the **dfoptim** package contains the Nelder-Mead routine written in R instead of C, it was easier and way faster to modify.

I added two things: 

- A data.frame containing all the updates down in the optimization procedure. So every single update can be accessed later on.
- A model argument containing the description whether to use the GEV/GP distribution (extreme value theory).

**Important note**: Since the later one is constituting the table containing the results, the package can only be used for functions requiring two ("gpd") or three ("gev") input parameters.
