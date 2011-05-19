class HermitianCode():
    def __init__(self, m, a=None):
        '''One-point AG code on Hermitian curve H_m with D=a*Q'''
        if not is_prime_power(m):
            raise Exception('m must be a prime power')
        self.m = m
        self.g = m*(m-1)/2
        field_size = m**2
        self.F = GF(field_size,'w')
        x, y, z = PolynomialRing(self.F, 3, 'xyz').gens()
        self.x, self.y, self.z = x, y, z
        f = x**(m+1) + y**m * z + y * z**m
        self.C = Curve(f)
        if a is None:
            self.a = 2*self.g
        else:
            self.a = a
        # t - number of errors that can be corrected with SV algorithm
        self.t = (self.a-3*self.g+1)/2
        self.L_D = self.L(self.a)
        self.L_A = self.L(self.t + self.g)
        self.L_B = self.L(self.t + 3*self.g + 1)
        self.L_C = self.L(self.a - self.t - self.g)
        self.S = Matrix(self.L_A).transpose() * Matrix(self.L_C)
        print self.S
        self.H()
        self.G()
        self.n = len(self.P)
        self.k = self.G().nrows()

    def __repr__(self):
        return 'AG code on Hermitian curve H_%d <%s> with D=%d*Q' % (self.m, self.C, self.a)

    def L(self, a):
        '''L(a*Q)'''
        if a <= 2*self.g - 2:
            raise Exception('''The degree of Q is not greater than 2*g - 2
                %d <= 2*%d - 2''' % (a, self.g))
        L_D_basis = []
        m = self.m
        x, y, z = self.x, self.y, self.z
        for i in range(0, m + 1):
            for j in range(0, a/(m+1) + 1):
                if i*m + j*(m+1) <= a:
                    #print 'x^%d y^%d / z^%d' % (i, j, i+j)
                    L_D_basis.append((x**i * y**j)/z**(i+j))
        if len(L_D_basis) != a + 1 - self.g:
            raise Exception('''The number of functions found does not
            satisfy the Riemann-Roch theorem
            Found: %d, needed %d''' % (len(L_D_basis), a + 1 - self.g))
        return L_D_basis

    def rational_points(self):
        try:
            return self._points
        except:
            self._points = self.C.rational_points()
            return self._points

    def H(self):
        try:
            return self._H
        except:
            points = self.rational_points()
            self.Q = points[1]
            self.P = [points[0]] + points[2:]
            rows = []
            for f in self.L_D:
                row = []
                for point in self.P:
                    row.append(f(point[0], point[1], point[2]))
                rows.append(row)
            # parity-check matrix
            self._H = Matrix(self.C.base_ring(), rows)
            return self._H

    def G(self):
        try:
            return self._G
        except:
            self._G = Matrix(self.C.base_ring(), self.H().transpose().kernel().basis())
            return self._G

    def encode(self, w):
        return vector(self.C.base_ring(), w) * self.G()

    def _apply(self, f, p):
        return f(p[0], p[1], p[2])

    def multiply(self, v, f):
        '''vector-function multiplication'''
        accumulator = 0
        for p_i, v_i in zip(self.P, v):
            accumulator += self._apply(f, p_i)*v_i
        return accumulator

    def decode(self, v):
        new_rows = []
        for row in self.S:
            new_row = []
            for f in row:
                new_row.append(self.multiply(v, f))
            new_rows.append(new_row)
        S = Matrix(self.C.base_ring(), new_rows)
        #print S
        theta = vector(self.L_A)*vector(S.kernel().basis()[0])
        error_positions = []
        for i, p in enumerate(self.P):
            if self._apply(theta, p) == 0:
                error_positions.append(i)
        error_value_system = []
        for f in self.L_B:
            row = []
            for i in error_positions:
                row.append(self._apply(f, self.P[i]))
            row.append(self.multiply(v, f))
            error_value_system.append(row)
        error_value_system = Matrix(self.C.base_ring(), error_value_system).echelon_form()
        #print error_value_system
        print 'Error position and value:'
        for i, pos in enumerate(error_positions):
            print pos, error_value_system[i][-1]
            v[pos] -= error_value_system[i][-1]
        decode_g = Matrix(self.C.base_ring(), self.G().rows() + [v]).transpose().echelon_form()
        decoded_message = [decode_g[i][-1] for i in range(0, decode_g.ncols()-1)]
        return decoded_message

def test(AG):
    print AG
    print AG.G()
    print 'We may add up to %d errors for SV algorithm' % AG.t
    for j in range(0, 10):
        w = list(AG.C.base_ring())[1]
        message = [w**(i+j) for i in range(0, AG.k)]
        print 'Original message:\n%s' % message
        v=AG.encode(message)
        print 'Encoded message:\n%s' % v
        for i in range(0, AG.t):
            v[i*2] += w**(i+j)
        print 'Message with errors:\n%s' % v
        print 'Decoding...'
        decoded = AG.decode(v)
        print 'Decoded message:\n%s' % decoded
        if message == decoded:
            print 'Decoded correctly'
        else:
            print 'Failed to decode correctly'
        #if raw_input() == 'q':
        #    break

if __name__ == '__main__':
    test(HermitianCode(2, 6))
    #test(HermitianCode(4, 40))

