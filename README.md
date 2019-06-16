# structure-factors
#calcualte structure factors (Sq) for polymer bulk system

该脚本用于计算聚合物本体自组装体系中不同相结构的结构参数

其中:
sq-calculate.py为单线程版本，q的值按照bin_size增加并求平均，比较粗糙，但是曲线也比较平滑；

sq-calculate-multiple-.py为多线程版本，其中：

sq-calculate-multiple-1.py的q值按照bin_size增加并求平均，比较粗糙，但是曲线也比较平滑；

sq-calculate-multiple-2.py的q值按照每一个bin_size中Sq值最高那个q作为这一个bin的q值，对LAM相效果比较好，而对HEX相一般；

sq-calculate-multiple-3.py的q值没有取平均，而是全部输入，使用时可以根据具体情况再选择前面两种。

input file:
    strfact.xyz
    md.gro
    
output file:
    sq.dat
    并且会在屏幕输出第一个主峰的q及对应的d=2pi/q
    
usage:
    python < sq-calculate-.py
    
