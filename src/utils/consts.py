c = 299792.4580   # [km/s] speed of light

lines = {
         'Halpha': (6562.808, -15.5, 15.5),   # Kuerster et al. (2003, A&A, 403, 1077)
         'Halpha': (6562.808, -40., 40.),
         'Haleft': (6562.808, -300., -100.),
         'Harigh': (6562.808, 100, 300),
         'Haleft': (6562.808, -500., -300.),
         'Harigh': (6562.808, 300, 500),
         'CaI':    (6572.795, -15.5, 15.5),   # Kuerster et al. (2003, A&A, 403, 1077)
         'CaH':    (3968.470, -1.09/2./3968.470*c, 1.09/2./3968.470*c), # Lovis et al. (2011, arXiv1107.5325)
         'CaK':    (3933.664, -1.09/2./3933.664*c, 1.09/2./3933.664*c), # Lovis et al.
         'CaIRT1': (8498.02, -15., 15.),      # NIST + my definition
         'CaIRT1a': (8492, -40, 40),          # My definition
         'CaIRT1b': (8504, -40, 40),          # My definition
         'CaIRT2': (8542.09, -15., 15.),      # NIST + my definition
         'CaIRT2a': (8542.09, -300, -200),    # My definition
         'CaIRT2b': (8542.09, 250, 350),      # My definition, +50 due telluric
         'CaIRT3': (8662.14, -15., 15.),      # NIST + my definition
         'CaIRT3a': (8662.14, -300, -200),    # NIST + my definition
         'CaIRT3b': (8662.14, 200, 300),      # NIST + my definition
         'NaD1':   (5889.950943, -15., 15.),  # NIST + my definition
         'NaD2':   (5895.924237, -15., 15.),  # NIST + my definition
         'NaDref1': (5885, -40, 40),          # My definition
         'NaDref2': ((5889.950+5895.924)/2, -40, 40),   # My definition
         'NaDref3': (5905, -40, 40)           # My definition
}
