import ee

global_ee = ee 

def calc_spek_rad(image):
    """1.Landsat 8 tasviri uchun spektral radiance hisoblash funksiyasi"""

    for x in range(1, 12):  # B1 dan B11 gacha  0 dan boshlanmaydi !
        qcal = image.select(f'B{x}')
        m_l = ee.Number(image.get(f'RADIANCE_MULT_BAND_{x}'))
        add = ee.Number(image.get(f'RADIANCE_ADD_BAND_{x}'))

        # Spektral radiance hisoblash
        radiance = qcal.multiply(m_l).add(add)

        # Yangi band qo‘shish
        image = image.addBands(radiance.rename(f'SR_{x}'))

    return image

def calc_ref(image):
        """2.Landsat 8 tasviri uchun Reflectivityni hisoblash funksiyasi"""
        # pi qiymati
        pi = ee.Number(3.141592)

        # dr ni hisoblash
        doy = ee.Number(ee.Date(image.get('system:time_start')).getRelative('day', 'year')).add(1)
        dr = ee.Number.expression(
            '1 + 0.033 * cos(DOY * 2 * PI / 365)',
            {
                'DOY': doy,
                'PI': pi
            }
        )

        # cos() ni hisoblash
        sun_elev = ee.Number(image.get('SUN_ELEVATION'))
        alfa = ee.Number(90).subtract(sun_elev)
        cos_alfa = alfa.multiply(pi).divide(180).cos()

        # esun qiymatlari
        esun_org = [1895.33, 2004.57, 1820.75, 1549.49, 951.76, 247.55, 85.46, 1723.88, 366.97, 1, 1]  # list 0 dan boshlanadi , dict kabi yozsa ham bo`ladi{1:1890, 2:2004, }

        # LOOP yozish
        for n in range(1, 12):  # 1 dan 11 gacha bandlar
            esun = ee.Number(esun_org[n-1])  # esun ni olish natija 0 sababi list 0 dan boshlandi
            band = image.select(f'SR_{n}')

            # Reflectivity ni hisoblash
            ref = image.expression(
                '(pi * band) / (esun * dr * cos_alfa)',
                {
                    'pi': pi,
                    'band': band,
                    'esun': esun,
                    'dr': dr,
                    'cos_alfa': cos_alfa
                }
            )

            #   bandni Nomi bilan qo‘shish
            image = image.addBands(ref.rename(f'REF_{n}'))

        return image
    # yakuniy natija REF_1, REF_2 ... kabi boladi
        
def calc_alb_toa(image):
    """3.Landsat 8 tasviri uchun Albedo at the Top Atmosphere (TOA) ni hisoblash funksiyasi"""

    # Landsat 8 uchun ESUN qiymatlari (B1-B11)
    esun_org = [1895.33, 2004.57, 1820.75, 1549.49, 951.76, 247.55, 85.46, 1723.88, 366.97, 1, 1]

    # jami ESUN ni hisoblash
    sum_esun = sum(esun_org[:10])  # Faqat B1-B10 ni hisobga olamiz, B11 = 1

    # Albedo TOA hisoblash uchun boshlang'ich qiymatlar
    weighted_sum = ee.Image(0)

    # Barcha bandlar bo‘yicha sikl
    for n in range(1, 11):  # B1-B10 uchun ishlaydi, B11 termal band (kerak emas)
        esun = ee.Number(esun_org[n - 1])  # `n - 1` qilish kerak, chunki Python indekslash 0 dan boshlanadi
        w = esun.divide(sum_esun)  # W_i = ESUN_i / ∑ESUN

        # REF bandni olish
        ref_band = image.select(f'REF_{n}')

        # TOA Albedo uchun hisoblash
        weighted_sum = weighted_sum.add(ref_band.multiply(w))

        # Yakuniy natija
        albedo_toa = weighted_sum.rename('ALBEDO_TOA')

    return image.addBands(albedo_toa)

def calc_albedo(image, srtm):
        """4.Landsat 8 tasviri uchun Albedo ni hisoblash funksiyasi"""
        tsw = ee.Image(0.75).add(ee.Image(2e-5).multiply(srtm))
        tsw = tsw.pow(2)
        path_radian = ee.Image.constant(0.03)  # Bastiaansen taklif qilgan qiymat
        toa = image.select('ALBEDO_TOA')
        albedo = (toa.subtract(path_radian)).divide(tsw)

        return image.addBands(albedo.rename('ALBEDO')) 

def calc_inshort_radiation(image, srtm):
        """5.Landsat 8 tasviri uchun Incoming Shortwave Radiationni hisoblash funksiyasi"""
        pi = ee.Number(3.141592)

        # tsw hisoblash
        tsw = ee.Image(0.75).add(ee.Image(2e-5).multiply(srtm))
        tsw = tsw.pow(2)

        # dr - Earth-Sun distance correction factor (IMAGE ga o‘giramiz)
        doy = ee.Number(ee.Date(image.get('system:time_start')).getRelative('day', 'year')).add(1)
        dr = ee.Image.constant(ee.Number(1).add(ee.Number(0.033).multiply(doy.multiply(ee.Number(2)).multiply(pi).divide(365).cos())))

        # cos() hisoblash
        sun_elev = ee.Number(image.get('SUN_ELEVATION'))
        alfa = ee.Number(90).subtract(sun_elev)
        cos_alfa = ee.Image.constant(alfa.multiply(pi).divide(180).cos())

        # G (quyosh doimiysi) IMAGE sifatida
        g = ee.Image.constant(1367)

        # Yakuniy natija IMAGE bo‘lishi kerak
        inshort_radiation = g.multiply(dr).multiply(cos_alfa).multiply(tsw)

        return image.addBands(inshort_radiation.rename('IN_SHORT_RAD'))
    
def calc_surface_ems(image):
        """6.Landsat 8 tasviri uchun Surface Emissivityni hisoblash funksiyasi"""

        nir = image.select('REF_5')
        red = image.select('REF_4')

        # NDVI hisoblash
        ndvi = image.normalizedDifference(['REF_5', 'REF_4']).rename('NDVI')

        # SAVI hisoblash
        l = ee.Image.constant(0.5)
        savi = image.expression(
            '1.5 * (nir - red) / (nir + red + l)',
            {
                'nir': nir,
                'red': red,
                'l': l
            }
        ).rename('SAVI')

        # LAI hisoblash
        lai = ee.Image.constant(0.69).subtract(savi).divide(0.59).log().divide(ee.Image.constant(0.91)).rename('LAI')

        # Surface emissivity E_NB va E_O
        e_nb = ee.Image(0.98).where(ndvi.gt(0).And(lai.lt(3)), ee.Image(0.97).add(ee.Image.constant(0.0033).multiply(lai))).rename('E_NB')
        e_o = ee.Image(0.98).where(ndvi.gt(0).And(lai.lt(3)), ee.Image(0.95).add(ee.Image.constant(0.01).multiply(lai))).rename('E_O')

        # Natijalarni tasvirga qo'shish
        return image.addBands([e_nb, e_o, lai, ndvi, savi])
    
def calc_surface_temp(image):
        """7.Landsat 8 tasviri uchun Surface Temperature hisoblash funksiyasi"""

        # K1 va K2 konstantalarni olish
        k_1 = ee.Number(image.get('K1_CONSTANT_BAND_10'))
        k_2 = ee.Number(image.get('K2_CONSTANT_BAND_10'))

        # Temperature bandini olish
        radiance = image.select('SR_10')

        # Surface Temperature (T_S) hisoblash
        t_s = image.expression(
            'K2 / log(ENB * K1 / Rc + 1)',
            {
                'K2': k_2,
                'ENB': image.select('E_NB'),
                'K1': k_1,
                'Rc': radiance
            }
        ).rename('T_S')

        # Natijani image ga qo'shish
        return image.addBands(t_s)

def calc_out_rad(image):
        """8.Landsat 8 tasviri uchun OUTGOING LONGWAVE RADIATION (OLR) hisoblash funksiyasi"""

        # Kerakli bandlarni tanlash
        t_s = image.select('T_S')  # Surface Temperature (Kelvin)
        e_o = image.select('E_O')  # Surface emissivity

        # Stefan-Boltzmann doimiysi (W/m^2K^4)
        stef_bolt = ee.Image.constant(5.67e-8)

        # Outgoing Longwave Radiation (OLR) hisoblash
        out_rad = stef_bolt.multiply(e_o).multiply(t_s.pow(4)).rename('OUT_RAD')

        # Natijani image ga qo‘shish
        return image.addBands(out_rad)
    
def calc_in_long(image, air_temp, srtm):
        """9.Landsat 8 tasviri uchun Incoming LONGWAVE RADIATION (OLR) hisoblash funksiyasi"""

        # tsw (atmosferik transmittans) ni hisoblash
        tsw = ee.Image(0.75).add(ee.Image(2e-5).multiply(srtm))
        a = (tsw.log().multiply(-1)).pow(0.09)

        # e_a (atmosferik emissivlik) ni hisoblash
        e_a = ee.Image.constant(0.85).multiply(a)

        # Stefan-Boltzmann doimiysi
        step_boltz = ee.Image.constant(5.67e-8)

        # Havo harorati (air_temp ni tanlash)
        t_air = air_temp.select('temperature_2m')
        t_air = t_air.pow(4)

        # Incoming Longwave Radiation (OLR) hisoblash
        in_long = step_boltz.multiply(e_a).multiply(t_air).rename('IN_LONG')

        # Natijani image ga qo'shish
        return image.addBands(in_long)

def calc_net_rad(image):
        """10.Landsat 8 tasviri uchun NET RADIATION (Rn) hisoblash funksiyasi"""
        albedo = image.select('ALBEDO')
        inshort_rad = image.select('IN_SHORT_RAD')
        out_rad = image.select('OUT_RAD')
        in_long = image.select('IN_LONG')
        e_o = image.select('E_O')

        # Rn ni hisoblash
        net_rad = image.expression(
            '(1-a)*RS - RL_T + RL_P * E_O',
            {
                'a': albedo,
                'RS': inshort_rad,
                'RL_T': out_rad,
                'RL_P': in_long,
                'E_O': e_o
            }
        ).rename('Rn')

        return image.addBands(net_rad)
    

