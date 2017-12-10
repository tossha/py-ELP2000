import math

# precession constant in J2000
_p = 5029.0966
_files = False

######################################################
#################  PUBLIC FUNCTIONS  #################
######################################################

#################  readFiles(path)
# Loads koefficients from files to memory
#
# Parameters:
#   path (str) — directory where the files are located, including trailing slash
def readFiles(path):
	global _files
	_files = []
	for i in range(36):
		_files.append([])
		with open(path + 'ELP' + str(i+1)) as file:
			isFirstLine = True
			for line in file:
				if isFirstLine:
					isFirstLine = False
					continue

				koeffs = {
					'i1': int(line[0:3]),
					'i2': int(line[3:6]),
					'i3': int(line[6:9]),
					'i4': int(line[9:12])
				}

				if i < 3:
					koeffs['A']  = float(line[14:27])
				elif i < 9:
					koeffs['i5'] = int  (line[12:15])
					koeffs['ph'] = float(line[16:25])
					koeffs['A']  = float(line[26:35])
				elif i < 21:
					koeffs['i5'] = int  (line[12:15])
					koeffs['i6'] = int  (line[15:18])
					koeffs['i7'] = int  (line[18:21])
					koeffs['i8'] = int  (line[21:24])
					koeffs['i9'] = int  (line[24:27])
					koeffs['i10']= int  (line[27:30])
					koeffs['i11']= int  (line[30:33])
					koeffs['ph'] = float(line[34:43])
					koeffs['A']  = float(line[44:53])
				else:
					koeffs['i5'] = int  (line[12:15])
					koeffs['ph'] = float(line[16:25])
					koeffs['A']  = float(line[26:35])

				_files[i].append(koeffs)

#################  getState(et)
# Calulates cartesian state vector (position and velocity) at given ephemeris time
#
# Parameters:
#   et (float) — ephemeris time (seconds past J2000 epoch)
# Returns:
#   (x, y, z, vx, vy, vz)
def getState(et):
	global _files
	if _files == False:
		print('Please load the files first: elp2000.readFiles(path)')
		return None

	model = _calcAll(et / 86400 / 36525)
	sinLon = math.sin(model['lon'])
	sinLat = math.sin(model['lat'])
	cosLon = math.cos(model['lon'])
	cosLat = math.cos(model['lat'])

	return (
		model['r'] * cosLon * cosLat,
		model['r'] * sinLon * cosLat,
		model['r'] * sinLat,
		(model['dr'] * cosLon * cosLat - model['r'] * sinLon * cosLat * model['dlon'] - model['r'] * cosLon * sinLat * model['dlat']) / 36525 / 86400,
		(model['dr'] * sinLon * cosLat + model['r'] * cosLon * cosLat * model['dlon'] - model['r'] * sinLon * sinLat * model['dlat']) / 36525 / 86400,
		(model['dr'] * sinLat + model['r'] * cosLat * model['dlat']) / 36525 / 86400
	)

######################################################
#################  PRIVATE FUNCTIONS  ################
######################################################

def _deg2rad(deg):
	return deg / 180 * math.pi

def _calc1_3(fileIdx, t):
	global _files
	D  = 297 * 3600 + 51 * 60 +  0.73512 + 1602961601.4603 * t -  5.8681 * t**2 + 0.006595 * t**3 - 0.00003184 * t**4
	dD = 1602961601.4603 - 5.8681 * 2 * t + 0.006595 * 3 * t**2 - 0.00003184 * 4 * t**3
	l_ = 357 * 3600 + 31 * 60 + 44.79306 +  129596581.0474 * t -  0.5529 * t**2 + 0.000147 * t**3
	dl_= 129596581.0474 - 0.5529 * 2 * t + 0.000147 * 3 * t**2
	l  = 134 * 3600 + 57 * 60 + 48.28096 + 1717915923.4728 * t + 32.3893 * t**2 + 0.051651 * t**3 - 0.00024471 * t**4
	dl = 1717915923.4728 + 32.3893 * 2 * t + 0.051651 * 3 * t**2 - 0.00024471 * 4 * t**3
	F  =  93 * 3600 + 16 * 60 + 19.55755 + 1739527263.0983 * t - 12.2505 * t**2 - 0.001021 * t**3 + 0.00000417 * t**4
	dF = 1739527263.0983 - 12.2505 * 2 * t - 0.001021 * 3 * t**2 + 0.00000417 * 4 * t**3

	value = 0
	derivative = 0

	for k in _files[fileIdx]:
		arg = _deg2rad((k['i1'] * D + k['i2'] * l_ + k['i3'] * l + k['i4'] * F) / 3600)
		dArg = k['i1'] * dD + k['i2'] * dl_ + k['i3'] * dl + k['i4'] * dF
		dArg = k['i1'] * dD + k['i2'] * dl_ + k['i3'] * dl + k['i4'] * dF
		sin = math.sin(arg)
		cos = math.cos(arg)
		value      += k['A'] * ( cos if fileIdx == 2 else sin)
		derivative += k['A'] * (-sin if fileIdx == 2 else cos) * _deg2rad(dArg/3600)

	return [value, derivative]

def _calc4_9(fileIdx, t):
	global _files
	D  = 297 * 3600 + 51 * 60 +  0.73512 + 1602961601.4603  * t
	dD = 1602961601.4603
	l_ = 357 * 3600 + 31 * 60 + 44.79306 +  129596581.0474  * t
	dl_= 129596581.0474
	l  = 134 * 3600 + 57 * 60 + 48.28096 + 1717915923.4728  * t
	dl = 1717915923.4728
	F  =  93 * 3600 + 16 * 60 + 19.55755 + 1739527263.0983  * t
	dF = 1739527263.0983
	W1 = 218 * 3600 + 18 * 60 + 59.95571 + 1732559343.73604 * t -  5.8883 * t**2 + 0.006604 * t**3 - 0.00003169 * t**4
	dW1= 1732559343.73604 - 5.8883 * 2 * t + 0.006604 * 3 * t**2 - 0.00003169 * 4 * t**3
	z = W1 + _p * t
	dz = dW1 + _p

	value = 0
	derivative = 0

	for k in _files[fileIdx]:
		arg = (k['i1'] * z + k['i2'] * D + k['i3'] * l_ + k['i4'] * l + k['i5'] * F) / 3600 + k['ph']
		dArg = k['i1'] * dz + k['i2'] * dD + k['i3'] * dl_ + k['i4'] * dl + k['i5'] * dF
		value      += k['A'] * math.sin(_deg2rad(arg))
		derivative += k['A'] * math.cos(_deg2rad(arg)) * _deg2rad(dArg/3600)

	return [value, derivative]

def _calc10_15(fileIdx, t):
	global _files
	# T  = 100 * 3600 + 27 * 60 + 59.22059 +  129597742.2758 * t -  0.0202 * t**2 + 0.000009 * t**3 + 0.00000015 * t**4
	D  = 297 * 3600 + 51 * 60 +  0.73512 + 1602961601.4603 * t
	dD = 1602961601.4603
	l  = 134 * 3600 + 57 * 60 + 48.28096 + 1717915923.4728 * t
	dl = 1717915923.4728
	F  =  93 * 3600 + 16 * 60 + 19.55755 + 1739527263.0983 * t
	dF = 1739527263.0983
	T  = 100 * 3600 + 27 * 60 + 59.22059 +  129597742.2758 * t
	dT =  129597742.2758
	Me = 252 * 3600 + 15 * 60 +  3.25986 + 538101628.68898 * t
	dMe= 538101628.68898
	V  = 181 * 3600 + 58 * 60 + 47.28305 + 210664136.43355 * t
	dV = 210664136.43355
	Ma = 355 * 3600 + 25 * 60 + 59.78866 +  68905077.59284 * t
	dMa=  68905077.59284
	J  =  34 * 3600 + 21 * 60 +  5.34212 +  10925660.42861 * t
	dJ =  10925660.42861
	S  =  50 * 3600 +  4 * 60 + 38.89694 +   4399609.65932 * t
	dS =   4399609.65932
	U  = 314 * 3600 +  3 * 60 + 18.01841 +   1542481.19393 * t
	dU =   1542481.19393
	N  = 304 * 3600 + 20 + 60 + 55.19575 +    786550.32074 * t
	dN =    786550.32074

	value = 0
	derivative = 0

	for k in _files[fileIdx]:
		arg = (k['i1'] * Me + k['i2'] * V + k['i3'] * T + k['i4'] * Ma + k['i5'] * J + k['i6'] * S + k['i7'] * U + k['i8'] * N + k['i9'] * D + k['i10'] * l + k['i11'] * F) / 3600 + k['ph']
		dArg = k['i1'] * dMe + k['i2'] * dV + k['i3'] * dT + k['i4'] * dMa + k['i5'] * dJ + k['i6'] * dS + k['i7'] * dU + k['i8'] * dN + k['i9'] * dD + k['i10'] * dl + k['i11'] * dF
		value      += k['A'] * math.sin(_deg2rad(arg))
		derivative += k['A'] * math.cos(_deg2rad(arg)) * _deg2rad(dArg/3600)

	return [value, derivative]

def _calc16_21(fileIdx, t):
	global _files
	D  = 297 * 3600 + 51 * 60 +  0.73512 + 1602961601.4603 * t
	dD = 1602961601.4603
	l_ = 357 * 3600 + 31 * 60 + 44.79306 +  129596581.0474 * t
	dl_=  129596581.0474
	l  = 134 * 3600 + 57 * 60 + 48.28096 + 1717915923.4728 * t
	dl = 1717915923.4728
	F  =  93 * 3600 + 16 * 60 + 19.55755 + 1739527263.0983 * t
	dF = 1739527263.0983
	T  = 100 * 3600 + 27 * 60 + 59.22059 +  129597742.2758 * t
	dT =  129597742.2758
	Me = 252 * 3600 + 15 * 60 +  3.25986 + 538101628.68898 * t
	dMe= 538101628.68898
	V  = 181 * 3600 + 58 * 60 + 47.28305 + 210664136.43355 * t
	dV = 210664136.43355
	Ma = 355 * 3600 + 25 * 60 + 59.78866 +  68905077.59284 * t
	dMa=  68905077.59284
	J  =  34 * 3600 + 21 * 60 +  5.34212 +  10925660.42861 * t
	dJ =  10925660.42861
	S  =  50 * 3600 +  4 * 60 + 38.89694 +   4399609.65932 * t
	dS =   4399609.65932
	U  = 314 * 3600 +  3 * 60 + 18.01841 +   1542481.19393 * t
	dU =   1542481.19393
	N  = 304 * 3600 + 20 + 60 + 55.19575 +    786550.32074 * t
	dN =    786550.32074

	value = 0
	derivative = 0

	for k in _files[fileIdx]:
		arg = (k['i1'] * Me + k['i2'] * V + k['i3'] * T + k['i4'] * Ma + k['i5'] * J + k['i6'] * S + k['i7'] * U + k['i8'] * D + k['i9'] * l_ + k['i10'] * l + k['i11'] * F) / 3600 + k['ph']
		dArg = k['i1'] * dMe + k['i2'] * dV + k['i3'] * dT + k['i4'] * dMa + k['i5'] * dJ + k['i6'] * dS + k['i7'] * dU + k['i8'] * dD + k['i9'] * dl_ + k['i10'] * dl + k['i11'] * dF
		value      += k['A'] * math.sin(_deg2rad(arg))
		derivative += k['A'] * math.cos(_deg2rad(arg)) * _deg2rad(dArg/3600)

	return [value, derivative]

def _calc22_36(fileIdx, t):
	global _files
	D  = 297 * 3600 + 51 * 60 +  0.73512 + 1602961601.4603 * t
	dD = 1602961601.4603
	l_ = 357 * 3600 + 31 * 60 + 44.79306 +  129596581.0474 * t
	dl_=  129596581.0474
	l  = 134 * 3600 + 57 * 60 + 48.28096 + 1717915923.4728 * t
	dl = 1717915923.4728
	F  =  93 * 3600 + 16 * 60 + 19.55755 + 1739527263.0983 * t
	dF = 1739527263.0983

	value = 0
	derivative = 0

	for k in _files[fileIdx]:
		arg = (k['i2'] * D + k['i3'] * l_ + k['i4'] * l + k['i5'] * F) / 3600 + k['ph']
		dArg = k['i2'] * dD + k['i3'] * dl_ + k['i4'] * dl + k['i5'] * dF
		value      += k['A'] * math.sin(_deg2rad(arg))
		derivative += k['A'] * math.cos(_deg2rad(arg)) * _deg2rad(dArg/3600)

	return [value, derivative]

def _calcFile(fileIdx, t):
	if fileIdx < 3:
		return _calc1_3(fileIdx, t)
	elif fileIdx < 6:
		return _calc4_9(fileIdx, t)
	elif fileIdx < 9:
		res = _calc4_9(fileIdx, t)
		return [res[0] * t, res[1] * t + res[0]]
	elif fileIdx < 12:
		return _calc10_15(fileIdx, t)
	elif fileIdx < 15:
		res = _calc10_15(fileIdx, t)
		return [res[0] * t, res[1] * t + res[0]]
	elif fileIdx < 18:
		return _calc16_21(fileIdx, t)
	elif fileIdx < 21:
		res = _calc16_21(fileIdx, t)
		return [res[0] * t, res[1] * t + res[0]]
	elif fileIdx < 24:
		return _calc22_36(fileIdx, t)
	elif fileIdx < 27:
		res = _calc22_36(fileIdx, t)
		return [res[0] * t, res[1] * t + res[0]]
	elif fileIdx < 33:
		return _calc22_36(fileIdx, t)
	elif fileIdx < 36:
		res = _calc22_36(fileIdx, t)
		return [res[0] * t**2, res[1] * t**2 + res[0] * 2 * t]

def _calcAll(t):
	lon  = 0
	dlon = 0
	lat  = 0
	dlat = 0
	r    = 0
	dr   = 0

	W1  = 218 * 3600 + 18 * 60 + 59.95571 + 1732559343.73604 * t - 5.8883 * t**2 + 0.006604 * t**3 - 0.00003169 * t**4
	dW1 = 1732559343.73604 - 5.8883 * 2 * t + 0.006604 * 3 * t**2 - 0.00003169 * 4 * t**3

	for i in range(12):
		f1 = _calcFile(i * 3    , t)
		f2 = _calcFile(i * 3 + 1, t)
		f3 = _calcFile(i * 3 + 2, t)
		lon  += f1[0]
		dlon += f1[1]
		lat  += f2[0]
		dlat += f2[1]
		r    += f3[0]
		dr   += f3[1]

	return {
		'lon':  _deg2rad(( lon +  W1) / 3600),
		'dlon': _deg2rad((dlon + dW1) / 3600),
		'lat':  _deg2rad( lat / 3600),
		'dlat': _deg2rad(dlat / 3600),
		'r':    r,
		'dr':   dr
	}
