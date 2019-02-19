# Description
A simple matrix calculator
# Functions
In addition to the basic operations, such as `==`, `!=`, `+`, `+=`, `-`, `-=`, `*` and `*=`, it also has the following functions.
|function|method name|
|:-:|:-:|
|get the element|(row,column)**(number start from zero)**|
|Row exchange|rowInterchange|
|Determinant|det|
|convert to Reduced Row-Echlon Form|rref|
|get the inverse Matrix|inv|
|get the transpose matrix|transpose|
|get the adjugate matrix|adjugate|
|get the characteristic polynomial|poly|
# Document
`伴随矩阵的多种求法.pdf` is a paper which I use the third method to calc the adjugate matrix
# Disadvantage
The time complexity of method `poly` is $O(n^2·n!)$ and method `adjugate` use many memory copy.
