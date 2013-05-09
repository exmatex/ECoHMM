
#ifndef COMDSWIFT_SWIFT
#define COMDSWIFT_SWIFT

(float virial_stress)
COMDSWIFT_runSim(string argv_string)
"comdswift" "0.0"
[ "set <<virial_stress>> [ COMDSWIFT_runSimString <<argv_string>> ]" ];

#endif
