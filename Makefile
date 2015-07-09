all:
  g++ lm.cpp `pkg-config --cflags --libs opencv` -o lm
clean:
  rm lm
  
