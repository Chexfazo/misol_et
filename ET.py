import ee
global_ee = ee

def calc_et(image):
        """ 1. Latent Heat Flux (λET) hisoblash funksiyasi """
        bands_needed = ['Rn', 'G', 'H']

        # Tasvirda kerakli bandlar borligini tekshiramiz
        missing_bands = [band for band in bands_needed if not image.bandNames().contains(band)]

        if missing_bands:
            raise ValueError(f"Quyidagi bandlar yetishmayapti: {missing_bands}")

        r_n = image.select('Rn')
        g = image.select('G')
        h = image.select('H')

        # ET = Rn - G - H
        et = r_n.subtract(g).subtract(h).rename('λET')

        return image.addBands(et)
    
def calc_ins_et(image):
        """ 2. Instantaneous ET (ETinst) ni hisoblash funksiyasi """
        et = image.select('λET')
        time = ee.Image.constant(3600)
        l = ee.Image.constant(2.5e6)
        et = et.multiply(time).divide(l).rename('ETinst')
        return image.addBands(et)
    
def calc_et_fraction(image):
        """ 3. ET fraction ni hisoblash funksiyasi """
        et = image.select('ETinst')
        h = image.select('H')
        et_fraction = et.divide(et.add(h)).rename('ET_F')
        return image.addBands(et_fraction)
    
def calc_ref_et(image, ET_0):
        """ 4. Reference ET_Fraction (Ref_ET) fraction ni hisoblash funksiyasi """
        et = image.select('ETinst')
        et_0 = ee.Image.constant(ET_0)
        et_fraction = et.divide(et_0).rename('Ref_ET')
        return image.addBands(et_fraction)
    
def calc_actual_et(image, ET_0):
        """ 5. Actual ET ni hisoblash funksiyasi """
        et = image.select('ETinst')
        et_0 = ee.Image.constant(ET_0)
        act_et = et.multiply(et_0).rename('Act_ET')
        return image.addBands(act_et)

import ee
global_ee = ee

def calc_et(image):
        """ 1. Latent Heat Flux (λET) hisoblash funksiyasi """
        bands_needed = ['Rn', 'G', 'H']

        # Tasvirda kerakli bandlar borligini tekshiramiz
        missing_bands = [band for band in bands_needed if not image.bandNames().contains(band)]

        if missing_bands:
            raise ValueError(f"Quyidagi bandlar yetishmayapti: {missing_bands}")

        r_n = image.select('Rn')
        g = image.select('G')
        h = image.select('H')

        # ET = Rn - G - H
        et = r_n.subtract(g).subtract(h).rename('λET')

        return image.addBands(et)
    
def calc_ins_et(image):
        """ 2. Instantaneous ET (ETinst) ni hisoblash funksiyasi """
        et = image.select('λET')
        time = ee.Image.constant(3600)
        l = ee.Image.constant(2.5e6)
        et = et.multiply(time).divide(l).rename('ETinst')
        return image.addBands(et)
    
def calc_et_fraction(image):
        """ 3. ET fraction ni hisoblash funksiyasi """
        et = image.select('ETinst')
        h = image.select('H')
        et_fraction = et.divide(et.add(h)).rename('ET_F')
        return image.addBands(et_fraction)
    
def calc_ref_et(image, ET_0):
        """ 4. Reference ET_Fraction (Ref_ET) fraction ni hisoblash funksiyasi """
        et = image.select('ETinst')
        et_0 = ee.Image.constant(ET_0)
        et_fraction = et.divide(et_0).rename('Ref_ET')
        return image.addBands(et_fraction)
    
def calc_actual_et(image, ET_0):
        """ 5. Actual ET ni hisoblash funksiyasi """
        et = image.select('ETinst')
        et_0 = ee.Image.constant(ET_0)
        act_et = et.multiply(et_0).rename('Act_ET')
        return image.addBands(act_et)
