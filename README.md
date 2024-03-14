# Multi-Threaded-RT

Here I use OpenMP to use multiple threads to render a scene. The image is divided into "tiles" of pixels(based on methods outlined in "Fundamentals of Computer Graphics). 
First, the number and allocation of tiles within the image to be rendered are calculated based on the given tile and image dimensions, and then the algorithm proceeds to allocate tiles to each thread accordingly.
Threads dynamically pick up and render a tile based on which tile is next to be rendered. 

Using 4 threads, there is a speed-up of at least 1.5x compared to rendering scenes using a single thread.
