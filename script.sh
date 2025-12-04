#!/bin/bash
#
#PBS -N 2uxj
#PBS -l mem=120GB
#PBS -l walltime=19:00:00
#PBS -l nodes=1:ppn=100
#PBS -j oe  # 将标准输出和错误输出合并到一个文件中。

# 使用指定网络接口用于 MPI 通信（如果集群需要）
export FI_TCP_IFACE=em2

# Calculate the number of processors allocated to this run.
NPROCS=`wc -l < $PBS_NODEFILE`

# Calculate the number of nodes allocated.
NNODES=`uniq $PBS_NODEFILE | wc -l`

echo "Running on host $(hostname)"
echo "Job started at $(date)"
echo "Allocated nodes: $NNODES"
echo "Total processors: $NPROCS"

# 设置库路径
export LD_LIBRARY_PATH=../ccp4_lib/lib:$LD_LIBRARY_PATH

# Optional, 苹果笔记本，动态链接库
#export DYLD_LIBRARY_PATH="/usr/local/Cellar/gcc/14.2.0_1/lib/gcc/14:../ccp4_lib/lib:$DYLD_LIBRARY_PATH"

# Optional, 学院服务器，启用mpicxx
#source /public/software/profile.d/mpi_intelmpi-2021.3.0.sh

# 进入工作目录
cd "$PBS_O_WORKDIR" || exit 1

# 编译项目
echo "Building project using Makefile..."
make clean
make -j "$NPROCS"

# 运行可执行程序（默认生成名为 'program'）
echo "Running the program..."
mpirun -np "$NPROCS" -f "$PBS_NODEFILE" ./program

echo "Job ended at $(date)"

