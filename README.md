# Lattice-based Verifiable Mixnet

Code accompannying the paper "Verifiable Mix-Nets and Distributed Decryption for Voting from Lattice-Based Assumptions" by Diego F. Aranha, Carsten Baum, Kristian Gjøsteen,
Tjerand Silde.

Depedencies are the NFLLib and FLINT 2.7.1 libraries.

For building the code, run `make` inside the source directory. This will build the binaries for `commit`, `encrypt`, `pianex` and `piaex` and `shuffle` to test and benchmark different modules of the code.
Respectively, these implement the commitment scheme, the distributed BGV cryptosystem, the two zero-knowledge proofs and the shuffle itself.

WARNING: This is an academic proof of concept, and in particular has not received code review. This implementation is NOT ready for any type of production use.
