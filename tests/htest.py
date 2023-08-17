
import re

hashmatcher = r"/[^a-f,^0-9]*([a-f, 0-9]{8})[^a-f,^0-9]*" # 8character hex hash
hashmatcher = r"([a-f, 0-9]{8})" # 8character hex
hashhashmatcher = r"[^a-f,^0-9]*([a-f, 0-9]{8})[^a-f,^0-9]*" # 8character hex hash

hmc = re.compile(hashmatcher)

na = '___abcdef78zzz'
nb = 'abcdef78'
n0 = 'adlskjfasl;dkfj;12a45f78ioqr0293uomnvca;lkn'
n1 = '2023-08-15_346699e0_2Dsearching_BH_simulation.csv'
n2 = 'fig2RedoExTime_346699e0.png'

tests = [na, nb, n0, n1, n2]

for t in tests:
    result = hmc.findall(t)
    print(f'{t:50} res: {result}')
    #print(result.groups())
