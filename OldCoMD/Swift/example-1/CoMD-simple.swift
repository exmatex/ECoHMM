
#include <builtins.swift>
#include <io.swift>

#include <CoMDSwift.swift>

main
{
  string s = "-f data/8k.inp.gz";
  int N = 3;
  foreach i in [0:N-1]
  {
    float virial_stress = COMDSWIFT_runSim(s);
    printf("Swift: virial_stress: %e", virial_stress);
  }
}
