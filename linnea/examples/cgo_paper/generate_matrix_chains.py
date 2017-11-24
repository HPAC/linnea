
import random
import math
import numpy


import itertools
from linnea.algebra import expression as ae
from linnea.algebra.equations import Equations
from linnea.algebra.transformations import simplify
from linnea.algebra.properties import Property as properties
from linnea.frontend.export import export
from linnea.code_generation import experiments as cge

def window(seq, n=2):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(itertools.islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result

def matrix_size():
    return random.randrange(50, 2001, 50)

_n = 0

def generate_operand(rows, columns):
    global _n
    _n += 1
    if columns == rows:
        if rows == 1:
            return # scalar

        operand = ae.Matrix('M{}'.format(_n), size=(rows, columns))
        operand.set_property(properties.FULL_RANK)
        # include "no property"
        if random.random() > 0.25:
            operand.set_property(random.choices([properties.DIAGONAL, properties.LOWER_TRIANGULAR, properties.UPPER_TRIANGULAR, properties.SYMMETRIC, properties.SPD])[0])
        return operand
    elif columns == 1:
        return ae.Vector('v{}'.format(_n), size=(rows, columns))
    elif rows == 1:
        return ae.Vector('v{}'.format(_n), size=(rows, columns))
    else:
        operand = ae.Matrix('M{}'.format(_n), size=(rows, columns))
        operand.set_property(properties.FULL_RANK)
        return operand

def matrix_chain_generator():
    while True:
        # for _ in range(20):
        #     expression_strategy.example()
        # generate_matrix_chain()

        length = random.randrange(4, 10)
        sizes = []

        if  0.95 <= random.random():
            sizes.append(1)
        else:
            sizes.append(matrix_size())

        for i in range(length):
            rand = random.random()
            if  0.6 <= rand < 0.95:
                # square
                sizes.append(sizes[i])
            elif 0.95 <= rand:
                # vector
                sizes.append(1)
            else:
                # non-square
                sizes.append(matrix_size())

        # print(sizes)
        operands = []
        restart = False
        for rows, columns in window(sizes):
            if rows == columns:
                if rows == 1:
                    restart = True
                    break
                    # operands.append(None) # error
                op = random.choices([ae.Identity, ae.Transpose, ae.Inverse, ae.InverseTranspose], weights=[2, 1, 1, 1])[0]
                operands.append(op(generate_operand(rows, columns)))
            elif rows == 1:
                operands.append(ae.Transpose(generate_operand(columns, rows)))
            elif columns == 1:
                operands.append(generate_operand(rows, columns))
            else:
                op = random.choices([ae.Identity, ae.Transpose], weights=[3, 1])[0]
                if op == ae.Transpose:
                    operands.append(op(generate_operand(columns, rows)))
                else:
                    operands.append(op(generate_operand(rows, columns)))
        if restart:
            continue

        # print(operands)
        expr = ae.Times(*operands)
        expr = simplify(expr)
        # print(expr)

        # print(expr.size)
        # print((sizes[0], sizes[-1]))
        if expr.has_property(properties.MATRIX):
            lhs = ae.Matrix('X'.format(_n), size=(sizes[0], sizes[-1]))
        elif expr.has_property(properties.VECTOR):
            lhs = ae.Vector('x'.format(_n), size=(sizes[0], sizes[-1]))

        yield Equations(ae.Equal(lhs, expr))

def spam():
    while True:
        yield 1

if __name__ == "__main__":

    from linnea import config

    config.set_language(config.Language.Julia)
    config.set_data_type(config.JuliaDataType.Float64)
    config.init()

    # from linnea.derivation.graph.constructive import DerivationGraph
    from linnea.derivation.graph.matrix_chain_derivation import DerivationGraph


    from timeit import default_timer as timer
    algos = 100
    iters = 10
    time_measurements = []
    paths=[]
    titles=[]
    measurement=True
    for seed in (100, 200, 300, 400, 500):#200, 300, 400):
        random.seed(seed)
        gen = matrix_chain_generator()
        for i in range(20):
            #if seed == 100 and i == 9:
            #    continue
            #if seed == 300 and (i >= 4 and i < 10):
            #    continue
            print("Seed {0}, chain {1}".format(seed, i))
            equations = next(gen)
            print(equations)
            if measurement:
                times = []
                for j in range(10):
                    start = timer()
                    graph = DerivationGraph(equations)
                    trace = graph.derivation(100, 10, verbose=False)
                    end = timer()
                    times.append(end - start)
                    if j == 0:
                        generated = graph.write_output(code=True,
                               pseudocode=True,
                               # output_name=args.input.split(".")[0],
                               output_name="seed_{}_matrix_chain_{}".format(seed, i),
                               operand_generator=True,
                               max_algorithms=100,
                               graph=True)
                if generated:
                    paths.append("seed_{}_matrix_chain_{}".format(seed, i))
                time_measurements.append([numpy.mean(times), numpy.std(times), numpy.min(times), numpy.max(times)])
                titles.append("seed_{}_matrix_chain_{}".format(seed, i))
            else:
                graph = DerivationGraph(equations)
                trace = graph.derivation(10, 6, verbose=False)
                generated = graph.write_output(code=True,
                               pseudocode=True,
                               # output_name=args.input.split(".")[0],
                               output_name="seed_{}_matrix_chain_{}".format(seed, i),
                               operand_generator=True,
                               max_algorithms=100,
                               graph=True)
                if generated:
                    paths.append("seed_{}_matrix_chain_{}".format(seed, i))
    if measurement:
        print(time_measurements)
        import csv
        with open('time_measurements', 'w') as f:
            writer = csv.writer(f,delimiter="\t")
            writer.writerow(["algorithm","Time","StdDev","Min","Max"])
            for row_title, data_row in zip(numpy.array(titles), numpy.array(time_measurements)):
                 writer.writerow([row_title] + data_row.tolist())
        cge.generate_cmake_script(paths)
                
        #rows = numpy.array(titles, dtype='|S20')[:, numpy.newaxis]
        #numpy.savetxt("generation_time", numpy.array(time_measurements), header="algorithm\tTime\tStdDev\tMin\tMax", fmt="%s %.18e %.18e %.18e %.18e")
        #with open('generation_time', 'w') as f:
        #numpy.savetxt(f, numpy.hstack((rows, numpy.array(time_measurements))), delimiter='\t', fmt='%s',header="algorithm\tTime\tStdDev\tMin\tMax")
    else:
        cge.generate_cmake_script(paths)



