# [Π: A Unified Framework for Computational Verifiable Secret Sharing](https://eprint.iacr.org/2023/1669)

This repository contains the implementation of the Verifiable Secret Sharing schemes Π_P​ and Π_LA​ Built by the General Framework Π, as presented in the paper titled [Π: A Unified Framework for Computational Verifiable Secret Sharing](https://eprint.iacr.org/2023/1669). The paper is publicly available on the [IACR eprint](https://eprint.iacr.org/2023/1669) archive and was authored by Karim Baghery.

## Overview
Current implementation is done using SageMath and the code includes the implementaiton of following schemes: 
- Shamir Secret Sharing 
- Pedersen Verifiable Secret Sharing
- ABCP23 Verifiable Secret Sharing scheme presented by Atapoor, Baghery, Cozzo, Pedersen  
- Our new Verifiable Secret Sharing Schemes Π_P​ and Π_LA 

## Running the Code

**WARNING:** This is an academic proof-of-concept prototype, and in particular has not received careful code review. This implementation is NOT ready for production use.

### Option 1: SageMath Terminal
A simple way to run our code is to copy the source code from this repository into the SageMath terminal and push the Enter press! 

### Option 2: Run in Online SageMath Terminal
Alternatively, you can run the code directly in an online SageMath terminal by clicking the following link:
[Run Online](https://sagecell.sagemath.org/????)

Please note that the online terminal may have memory limitations and cannot support more than n=512 parties due to these constraints. 

## License

This library is licensed under either of the following licenses, at your discretion.

 * Apache License Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)
