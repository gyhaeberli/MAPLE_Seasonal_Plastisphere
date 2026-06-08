This is the for Raw data structure

General library structure:

4 EUK (18S) libraries:

- GAB_EUK1: 2 plates → nested (one samples name is 2 samples with same position).
    
    
    | Sample plates | P1 | P2 |
    | --- | --- | --- |
    | Adapterama | set1 (A1,B2,C3,D4) | set2 (E5, F6, G7, H8) |
    | iNEXT | set1 | set1 |
- GAB_EUK2: 2 plates → P4 has a mixed set, but iNext differentiable
    
    
    | Sample plates | P3 | P4 |
    | --- | --- | --- |
    | Adapterama | set1 (A1,B2,C3,D4) | mix set (A1, F6, G7, D4) |
    | iNEXT | set1 | set2 |
- GAB_EUK3: 2 plates → iNext difference!
    
    
    | Sample plates | P5 | P6 |
    | --- | --- | --- |
    | Adapterama | set1 (A1,B2,C3,D4) | set1 (A1,B2,C3,D4) |
    | iNEXT | set1 | set2 |
- GAB_EUK4: 1 plate
    
    
    | Sample plates | P7 |
    | --- | --- |
    | Adapterama | set1 (A1,B2,C3,D4) |
    | iNEXT | set1 |