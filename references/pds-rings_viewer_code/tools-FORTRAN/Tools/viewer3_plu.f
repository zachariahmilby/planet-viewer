c***********************************************************************
c***********************************************************************

        program VIEWER3_PLU

        implicit                none
        include                 'tools.inc'

        integer                 NPLANET, PLANET_ID, BARYCENTER_ID
        parameter               (NPLANET = 9)
        parameter               (PLANET_ID = 100 * NPLANET + 99)
        parameter               (BARYCENTER_ID = NPLANET)

        character*(*)           PLANET_NAME, STARLIST, LONDIR
        parameter               (PLANET_NAME = 'Pluto')
        parameter               (STARLIST = 'starlist_plu.txt')
        parameter               (LONDIR = 'W')

        double precision        RPLANET, PC_DISTANCE
        parameter               (RPLANET = 1153.d0)
        parameter               (PC_DISTANCE = 19571.d0)

        integer                 NEW_HORIZONS_ID
        parameter               (NEW_HORIZONS_ID = -98)

c Moon parameters
        integer                 NMOONS
        parameter               (NMOONS = 5)
        integer                 MOON_IDS(0:NMOONS)/
     &                                  999,
     &                                  901, 902, 903, 904, 905 /
        character*10            MOON_NAMES(0:NMOONS)/
     &                                  'Pluto',
     &                                  'Charon', 'Nix', 'Hydra',
     &                                  'Kerberos', 'Styx' /
        logical                 MOON_IS_IRREGULAR(0:NMOONS)
     &                                  / .FALSE., NMOONS*.FALSE. /

c Ring constants
        integer                 NRINGS
        parameter               (NRINGS = 5)

        double precision        BARYCENTER_OFFSET
        parameter               (BARYCENTER_OFFSET = 2035.d0)

        double precision        ring_rads(NRINGS)/ 19571.d0,
     &                                             48703.d0,
     &                                             64754.d0,
     &                                             57771.d0,
     &                                             42608.d0 /
        double precision        ring_elevs(NRINGS)/     NRINGS * 0.d0 /
        double precision        ring_eccs(NRINGS)/      NRINGS * 0.d0 /
        double precision        ring_incs(NRINGS)/      NRINGS * 0.d0 /
        double precision        ring_peris(NRINGS)/     NRINGS * 0.d0 /
        double precision        ring_nodes(NRINGS)/     NRINGS * 0.d0 /
        double precision        ring_offsets(3,NRINGS)/ NRINGS * 0.d0,
     &                                                  NRINGS * 0.d0,
     &                                                  NRINGS * 0.d0 /

        logical                 ring_flags(NRINGS)/ 5 * .FALSE. /
        logical                 ring_opaqs(NRINGS)/ 5 * .FALSE. /
        logical                 ring_dashed(NRINGS)/ 5 * .TRUE. /

        integer                 NARCS
        parameter               (NARCS = 0)
        integer                 arc_rings/ 0 /
        double precision        arc_minlons/ 0.d0 /, arc_maxlons/ 0.d0 /
        double precision        arc_widthpts/ 0.d0 /
        logical                 arc_flags/ .FALSE. /

c Planet-specific declarations
        double precision        bary_dpv(6)
        integer                 imodel

c***********************************************************************

c Parameter parsing variables
        character*256           string, string2, time_string
        character*40            keyword
        integer                 i, isc, obs_isc, ivalue, iparen, jparen,
     &                          dutc
        double precision        secs, radius, lat, lon, alt, angle,
     &                          sun_range, obs_range, value
        logical                 status, found, blank_disks, isright

c Geometry variables
        double precision        planet_ra, planet_dec, ra, dec, dra,
     &                          ddec, cosdec, phase, obs_b, sun_b,
     &                          sun_db, obs_long, sun_long, obs_dist,
     &                          sun_dist, obs_pv(6), planet_dpv(6), dt,
     &                          subobs_lon, subobs_lat,
     &                          subsol_lon, subsol_lat
        character*18            ra_string, dec_string, litstr
        logical                 isdark

c Function declarations
        integer                 LASTNB
        double precision        RPD, DPR, CLIGHT, PI,
     &                          FJUL_ETofTAI, FJUL_TAIofDUTC
        logical                 RSPK_LoadFiles, RSPK_LoadSC,
     &                          FJUL_ParseDT, WWW_GetKey, WWW_GetKeys,
     &                          WWW_GetKeyUnsanitized,
     &                          ParseAngle, ReadStars

c Useful constants
        double precision        MAXSECS, AU
        parameter               (MAXSECS = 360.d0 * 60.d0 * 60.d0)
        parameter               (AU = 149597870.7d0)

c General plot parameters
        double precision        obs_time, fov, center_ra, center_dec
        integer                 center_body_id
        character*80            center_body_name, title, filename

        integer                 ncaptions, MAXCAPTIONS
        parameter               (MAXCAPTIONS = 8)
        character*80            lcaptions(MAXCAPTIONS)
        character*80            rcaptions(MAXCAPTIONS)

        logical                 moon_flags(0:NMOONS)
        double precision        moon_labelpts, moon_diampts,
     &                          meridian_pts
        integer                 sc_trajectory
        integer                 ring_method/0/

c Star parameters
        integer                 MAX_NSTARS
        parameter               (MAX_NSTARS = 200)

        integer                 nstars
        double precision        star_ras(MAX_NSTARS+10),
     &                          star_decs(MAX_NSTARS+10)
        character*20            star_names(MAX_NSTARS+10)

        logical                 star_labels
        double precision        STAR_DIAMPTS
        parameter               (STAR_DIAMPTS = 24.d0)

c***********************************************************************
c Initialize
c***********************************************************************

c Reduce SPICE error messages
        call ERRPRT('SET', 'NONE, SHORT, LONG, EXPLAIN')

c Disable tracing
        call TRCOFF

c***********************************************************************
c Summarize request
c***********************************************************************

10      format(a)
        write(*,10) 'Input Parameters'
        write(*,10) '----------------'
        write(*,10) ' '

c Observation time
        status = WWW_GetKey('time', string)
        write(*,10) '  Observation time: ' // string(1:LASTNB(string))

c Ephemeris selection
        status = WWW_GetKey('ephem', string)
        write(*,10) '         Ephemeris: ' // string(5:LASTNB(string))

c Field of view
        status = WWW_GetKey('fov', string)
        status = WWW_GetKey('fov_unit', string2)
        if (string .eq. ' ' .and. index(string2,'FOV') .gt. 0)
     &          string = '1'
        write(*,10) '     Field of view: ' //
     &              string(1:LASTNB(string)) //
     &              ' (' // string2(1:LASTNB(string2)) // ')'

c Diagram center
        status = WWW_GetKey('center', string)
        if (string .eq. 'body') then
            status = WWW_GetKey('center_body', string)
            write(*,10) '    Diagram center: ' //
     &                  string(1:LASTNB(string))
        else if (string .eq. 'ansa') then
            status = WWW_GetKey('center_ansa', string)
            status = WWW_GetKey('center_ew',   string2)
            write(*,10) '    Diagram center: ' //
     &                  string(1:LASTNB(string)) // ' ' //
     &                  string2(1:LASTNB(string2)) // ' ansa'
        else if (string .eq. 'J2000') then
            status = WWW_GetKey('center_ra', string)
            status = WWW_GetKey('center_ra_type', string2)
            write(*,10) '    Diagram center: RA  = ' //
     &                  string(1:LASTNB(string)) // ' ' //
     &                  string2(1:LASTNB(string2))
            status = WWW_GetKey('center_dec', string)
            write(*,10) '                    Dec = ' //
     &                  string(1:LASTNB(string))
        else
            status = WWW_GetKey('center_star', string)
            write(*,10) '    Diagram center: Star = ' //
     &                  string(1:LASTNB(string))
        end if

c Viewpoint
        status = WWW_GetKey('viewpoint', string)
        if (string .eq. 'latlon') then
            status = WWW_GetKey('latitude',  string)
            write(*,10) '         Viewpoint: Lat = ' //
     &                  string(1:LASTNB(string))  // ' (deg)'
            status = WWW_GetKey('longitude', string)
            status = WWW_GetKey('lon_dir',   string2)
            write(*,10) '                    Lon = ' //
     &                  string(1:LASTNB(string))  // ' (deg ' //
     &                  string2(1:LASTNB(string2))  // ')'
            status = WWW_GetKey('altitude',  string)
            write(*,10) '                    Alt = ' //
     &                  string(1:LASTNB(string)) // ' (m)'
        else
            status = WWW_GetKey('observatory', string)
            status = WWW_GetKey('sc_trajectory', string2)
            if (status) then
                string(LASTNB(string)+1:) =
     &                  ' (' // string2(5:LASTNB(string2)) // ')'
            end if

            write(*,10) '         Viewpoint: ' //
     &                  string(1:LASTNB(string))
        end if

c Moon selection
        status = WWW_GetKey('moons', string)
        write(*,10) '    Moon selection: ' // string(5:LASTNB(string))

        status = WWW_GetKey('moremoons', string)
        if (status) then
            write(*,10) '                    ' //
     &                  string(1:LASTNB(string))
        end if

c Ring selection
        status = WWW_GetKey('rings', string)
        write(*,10) '    Ring selection: ' // string(1:LASTNB(string))

c Standard stars
        status = WWW_GetKey('standard', string)
        if (.not. status) string = 'No'
        write(*,10) '    Standard stars: ' // string(1:LASTNB(string))

c Additional star
        status = WWW_GetKey('additional', string)
        if (.not. status) then
            write(*,10) '   Additional star: No'
        else
            status = WWW_GetKey('extra_name', string)
            write(*,10) '   Additional star: ' //
     &                  string(1:LASTNB(string))
            status = WWW_GetKey('extra_ra', string)
            status = WWW_GetKey('extra_ra_type', string2)
            write(*,10) '                    RA  = ' //
     &                  string(1:LASTNB(string)) // ' ' //
     &                  string2(1:LASTNB(string2))
            status = WWW_GetKey('extra_dec', string)
            write(*,10) '                    Dec = ' //
     &                  string(1:LASTNB(string))
        end if

c Other bodies
        status = WWW_GetKeys('other', string)
        if (.not. status) then
            write(*,10) '      Other bodies: None'
        else
            write(*,10) '      Other bodies: ' //
     &                  string(1:LASTNB(string))
            do 50 i = 1, 9999
                status = WWW_GetKeys('other', string)
                if (.not. status) goto 51
                write(*,10) '                    ' //
     &                  string(1:LASTNB(string))
50          continue
51          continue
        end if

c Title
        status = WWW_GetKey('title', string)
        write(*,10) '             Title: "' // string(1:LASTNB(string))
     &                                      // '"'

c Moon labels
        status = WWW_GetKey('labels', string)
        write(*,10) '       Moon labels: ' // string(1:LASTNB(string))

c Moon enlargement
        status = WWW_GetKey('moonpts', string)
        write(*,10) '  Moon enlargement: ' // string(1:LASTNB(string))
     &                  // ' (points)'

c Blank disks
        status = WWW_GetKey('blank', string)
        write(*,10) '       Blank disks: ' // string(1:LASTNB(string))

c Highlight prime meridians
        status = WWW_GetKey('meridians', string)
        write(*,10) '   Prime meridians: ' // string(1:LASTNB(string))

        write(*,10) ' '

c***********************************************************************
c Get time
c***********************************************************************

c Initialize leap seconds table
        call FJUL_InitLeaps(JULIAN_LEAPSECS)

c Extract time value
        keyword = 'OBSERVATION_TIME'
        status = WWW_GetKey('time', string)
        time_string = string

c Parse and convert to UTC
        status = FJUL_ParseDT(string, ' ', dutc, secs)
        if (.not. status) goto 9999

c Convert to ephemeris time
        obs_time = FJUL_ETofTAI(FJUL_TAIofDUTC(dutc) + secs)

c***********************************************************************
c Check need for spacecraft ephemerides
c
c This must be done BEFORE calling RSPK_LoadFiles(), because otherwise
c a spacecraft ephemeris will take precedence over the user's selection.
c***********************************************************************

c Check the observatory
        keyword = 'VIEWPOINT'
        status = WWW_GetKey('observatory', string)
        obs_isc = 0
        do isc = 1, NSPACECRAFTS
            if (string .eq. SC_NAMES(isc)) then
                obs_isc = isc
            end if
        end do

        keyword = 'SPACECRAFT_TRAJECTORY'
        sc_trajectory = 0
        status = WWW_GetKey('sc_trajectory', string)
        if (status) then
            read(string(1:4),*,err=9999) sc_trajectory
        end if

c Check the other bodies
        do 60 i = 1, 9999
            status = WWW_GetKeys('other', string)
            if (.not. status) goto 61

            do isc = 1, NSPACECRAFTS
              if (string .eq. SC_NAMES(isc) .and.
     &            isc .ne. obs_isc) then
                    status = RSPK_LoadSC(1, SC_IDS(isc), NPLANET,
     &                                   0, .FALSE.)
                    if (.not. status) goto 9999
              end if
            end do
60      continue
61      continue

c Load spacecraft ephemeris if necessary
        if (obs_isc .ne. 0) then
            status = RSPK_LoadSC(1, SC_IDS(obs_isc), NPLANET,
     &                           sc_trajectory, .TRUE.)
            if (.not. status) goto 9999
        end if

c***********************************************************************
c Set up ephemeris
c***********************************************************************

c Load SPICE files
        keyword = 'EPHEMERIS'
        status = WWW_GetKey('ephem', string)
        read(string(1:4), *, err=9999) ivalue

        status = RSPK_LoadFiles(1, NPLANET, ivalue)
        if (.not. status) goto 9999

c***********************************************************************
c Set viewpoint
c***********************************************************************

        keyword = 'VIEWPOINT'
        status = WWW_GetKey('viewpoint', string)

        if (.not. status) then
            call RSPK_SetObsId(EARTH_ID)
        else if (string .eq. 'observatory') then
            keyword = 'OBSERVATORY'
            status = WWW_GetKey('observatory', string)
            if (.not. status) goto 9999

            found = .FALSE.
            do isc = 1, NSPACECRAFTS
                if (string .eq. SC_NAMES(isc)) then
                    call RSPK_SetObsId(SC_CODES(isc))
                    found = .TRUE.
                end if
            end do

            if (.not. found) then
                if (string(1:5) .eq. 'Earth') then
                    call RSPK_SetObsId(EARTH_ID)
                else
                    iparen = index(string, '(')
                    jparen = index(string, ')')
                    read(string(iparen+1:jparen-1), *, err=9999)
     &                                              lat, lon, alt
                    call RSPK_SetObsId(EARTH_ID)
                    call RSPK_SetObs(lat, lon, alt)
                end if
            end if
        else if (string .eq. 'latlon') then
            keyword = 'LATITUDE'
            status = WWW_GetKey('latitude', string)
            status = ParseAngle(string, lat)
            if (.not. status) goto 9999

            keyword = 'LONGITUDE'
            status = WWW_GetKey('longitude', string)
            status = ParseAngle(string, lon)
            if (.not. status) goto 9999

            keyword = 'ALTITUDE'
            status = WWW_GetKey('altitude', string)
            read(string,*,err=9999) alt

            keyword = 'LONGITUDE_EW'
            status = WWW_GetKey('lon_dir', string)
            if (string .eq. 'east') then
                continue
            else if (string .eq. 'west') then
                lon = -lon
            else
                goto 9999
            end if

            call RSPK_SetObsId(EARTH_ID)
            call RSPK_SetObs(lat, lon, alt)
        else
            goto 9999
        end if

c***********************************************************************
c Set the field of view
c
c This must be handled after the ephemeris time is established and after
c RSPK_LoadFiles has been called, because it might have to determine
c the distance to the planet.
c***********************************************************************

c Get the field of view
        keyword = 'FIELD_OF_VIEW'
        status = WWW_GetKey('fov', string)
        status = WWW_GetKey('fov_unit', string2)
        if (string .eq. ' ' .and. index(string2,'FOV') .gt. 0)
     &          string = '1'
        if (string .eq. ' ') goto 9999
        read(string, *, err=9999) fov
        if (fov .le. 0.) goto 9999

c Get the field of view units
        keyword = 'FIELD_OF_VIEW_UNIT'
        status = WWW_GetKey('fov_unit', string)

        if (string(1:3) .eq. 'sec') then
            fov = fov/3600.d0 * RPD()
        else if (string .eq. 'degrees') then
            fov = fov * RPD()
        else if (string(1:8) .eq. 'millirad' .or.
     &           string(1:4) .eq. 'mrad') then
            fov = fov * 1.d-3
        else if (string(1:8) .eq. 'microrad' .or.
     &           string(1:4) .eq. 'urad') then
            fov = fov * 1.d-6
        else if (index(string, 'radi') .gt. 0) then
            call RSPK_Ranges(obs_time, sun_range, obs_range)
            fov = fov * asin(RPLANET / obs_range)
        else if (string(1:9) .eq. 'kilometer' .or.
     &           string(1:2) .eq. 'km') then
            call RSPK_Ranges(obs_time, sun_range, obs_range)
            fov = fov * asin(1.d0 / obs_range)
        else if (index(string, 'Charon') .gt. 0) then
            call RSPK_Ranges(obs_time, sun_range, obs_range)
            fov = fov * asin(PC_DISTANCE / obs_range)
        else if (index(string, 'LORRI') .gt. 0) then
            fov = fov * 0.2907d0 * RPD()
        else
            goto 9999
        end if

        fov = min(fov, PI()/2.d0)

c***********************************************************************
c Locate center of diagram
c***********************************************************************

c Get the central body ID
        center_body_id = PLANET_ID
        center_body_name = PLANET_NAME

        keyword = 'DIAGRAM_CENTER'
        status = WWW_GetKey('center', string)

c***************************************
c Case #1: A body
c***************************************

        if (string .eq. 'body') then
            keyword = 'CENTER_BODY'
            status = WWW_GetKey('center_body', string)

            if (string .eq. 'Barycenter') then
                center_body_id = PLANET_ID
                center_body_name = PLANET_NAME
                call RSPK_BodyRaDec(obs_time, BARYCENTER_ID,
     &                              center_ra, center_dec)
            else
              if (string .eq. PLANET_NAME) then
                center_body_id = PLANET_ID
                center_body_name = PLANET_NAME
              else
                iparen = index(string,'(')
                jparen = index(string,')')
                if (iparen .gt. 0) then
                    read(string(iparen+2:jparen-1), *, end=9999) ivalue
                    if (ivalue .le. 0) goto 9999
                    center_body_id = ivalue + 100 * NPLANET
                    center_body_name = string(1:iparen-2)
                else
                    goto 9999
                end if
              end if

              call RSPK_BodyRaDec(obs_time, center_body_id,
     &                            center_ra, center_dec)
            end if

c***************************************
c Case #2: A ring ansa
c***************************************

c Check for a ring
        else if (string .eq. 'ansa') then
            keyword = 'CENTER_ANSA'
            status = WWW_GetKey('center_ansa', string)

            if (index(string,'Nix') .gt. 0) then
                radius = ring_rads(2)
            else if (index(string,'Hydra') .gt. 0) then
                radius = ring_rads(3)
            else if (index(string,'Kerberos') .gt. 0) then
                radius = ring_rads(4)
            else if (index(string,'Styx') .gt. 0) then
                radius = ring_rads(5)
            else
                goto 9999
            end if

            keyword = 'CENTER_EAST_WEST'
            status = WWW_GetKey('center_ew', string)

            if (string .eq. 'west') then
                isright = .TRUE.
            else if (string .eq. 'east') then
                isright = .FALSE.
            else
                goto 9999
            end if

            call RSPK_AnsaRaDec(obs_time, radius, isright,
     &                  center_ra, center_dec)

c***************************************
c Case #3: J2000 coordinates
c***************************************

c Check for a J2000 ra & dec pair
        else if (string .eq. 'J2000') then

c Read three RA values
            keyword = 'RIGHT_ASCENSION'
            status = WWW_GetKey('center_ra', string)
            status = ParseAngle(string, angle)
            if (.not. status) goto 9999

            status = WWW_GetKey('center_ra_type', string)
            if (string(1:1) .eq. 'd' .or. string(1:1) .eq. 'D')
     &          angle = angle / 15.d0

            center_ra = angle * 15.d0 * RPD()

c Read three Dec values
            keyword = 'DECLINATION'
            status = WWW_GetKey('center_dec', string)
            status = ParseAngle(string, angle)
            if (.not. status) goto 9999
            center_dec = angle * RPD()

c***************************************
c Case #4: A Named star
c***************************************

        else if (string .eq. 'star') then
            keyword = 'STANDARD_STARLIST'
            status = ReadStars(1, STARLIST_PATH // STARLIST,
     &                            nstars, MAX_NSTARS,
     &                            star_ras, star_decs, star_names)
            if (.not. status) goto 9999

            keyword = 'CENTER_STAR'
            status = WWW_GetKey('center_star', string)
            do i = 1, nstars
                if (star_names(i) .eq. string) then
                    center_ra = star_ras(i)
                    center_dec = star_decs(i)
                    goto 71
                end if
            end do
            goto 9999

71          continue

c Unrecognized option
        else
                goto 9999
        end if

c***********************************************************************
c Moon selection
c***********************************************************************

        keyword = 'MOON_SELECTION'
        status = WWW_GetKey('moons', string)
        if (status) then
            read(string(1:3),*,err=9999) ivalue
            if (ivalue .le. 0) goto 9999
        else
            ivalue = PLANET_ID - 1
        end if

        moon_flags(0) = .TRUE.
        do i = 1, NMOONS
            moon_flags(i) = (MOON_IDS(i) .le. ivalue .and.
     &                       .not. MOON_IS_IRREGULAR(i))
        end do

c Select Irregulars
        status = WWW_GetKey('moremoons', string)
        if (string .ne. ' ') then
            do i = 1, NMOONS
                if (MOON_IS_IRREGULAR(i)) moon_flags(i) = .TRUE.
            end do
        end if

c***********************************************************************
c Ring selection
c***********************************************************************

        keyword = 'RING_SELECTION'
        status = WWW_GetKey('rings', string)

        if (index(string,'Charon')   .gt. 0) ring_flags(1) = .TRUE.
        if (index(string,'Nix')      .gt. 0) ring_flags(2) = .TRUE.
        if (index(string,'Hydra')    .gt. 0) ring_flags(3) = .TRUE.
        if (index(string,'Kerberos') .gt. 0) ring_flags(4) = .TRUE.
        if (index(string,'Styx')     .gt. 0) ring_flags(5) = .TRUE.

c***********************************************************************
c Background objects
c***********************************************************************

c Standard stars
        status = WWW_GetKey('standard', string)
        if (string .eq. 'Yes') then
            keyword = 'STANDARD_STARLIST'
            status = ReadStars(1, STARLIST_PATH // STARLIST,
     &                            nstars, MAX_NSTARS,
     &                            star_ras, star_decs, star_names)
                if (.not. status) goto 9999
        else
            nstars = 0
        end if

c Additional star
        status = WWW_GetKey('additional', string)
        if (string .eq. 'Yes') then
            nstars = nstars + 1

            status = WWW_GetKey('extra_name', star_names(nstars))

            keyword = 'STAR_RA'
            status = WWW_GetKey('extra_ra', string)
            status = ParseAngle(string, angle)
            if (.not. status) goto 9999

            status = WWW_GetKey('extra_ra_type', string)
            if (string(1:1) .eq. 'd' .or. string(1:1) .eq. 'D')
     &          angle = angle / 15.d0
            star_ras(nstars) = angle * 15.d0 * RPD()

            keyword = 'STAR_DEC'
            status = WWW_GetKey('extra_dec', string)
            status = ParseAngle(string, angle)
            if (.not. status) goto 9999
            star_decs(nstars) = angle * RPD()
        end if

c Other bodies
        keyword = 'OTHER_BODIES'
        do 500 i = 1, 9999
            status = WWW_GetKeys('other', string)
            if (.not. status) goto 501

            nstars = nstars + 1
            if (string .eq. 'Sun') then
                call RSPK_BodyRaDec(obs_time, SUN_ID,
     &                  star_ras(nstars), star_decs(nstars))
                star_names(nstars) = 'Sun'
            else if (string .eq. 'Anti-Sun') then
                call RSPK_AntiSun(obs_time, PLANET_ID,
     &                  star_ras(nstars), star_decs(nstars))
                star_names(nstars) = 'Anti-Sun'
            else if (string .eq. 'Earth') then
                call RSPK_BodyRaDec(obs_time, EARTH_ID,
     &                  star_ras(nstars), star_decs(nstars))
                star_names(nstars) = 'Earth'
            else if (string .eq. 'Barycenter') then
                call RSPK_BodyRaDec(obs_time, BARYCENTER_ID,
     &                  star_ras(nstars), star_decs(nstars))
                star_names(nstars) = string
            else
                found = .FALSE.
                do isc = 1, NSPACECRAFTS
                    if (string .eq. SC_NAMES(isc)) then
                        call RSPK_BodyRaDec(obs_time, SC_CODES(isc),
     &                      star_ras(nstars), star_decs(nstars))
                        star_names(nstars) = string
                        found = .TRUE.
                    end if
                end do

                if (.not. found) goto 9999
            end if
500     continue
501     continue

c***********************************************************************
c Get size of labels
c***********************************************************************

        keyword = 'MOON_LABELS'
        status = WWW_GetKey('labels', string)

        if (status) then
            iparen = index(string, '(')
            jparen = index(string, ' points)')
            read(string(iparen+1:jparen-1), *, err=9999) moon_labelpts
            if (moon_labelpts .le. 0.) goto 9999
        else
            moon_labelpts = 0.d0
        end if

c Label stars only if we're also labeling moons
        star_labels = (moon_labelpts .gt. 0.d0)

c***********************************************************************
c Get minimum disk size
c***********************************************************************

        keyword = 'MOON_ENLARGEMENT'
        status = WWW_GetKey('moonpts', string)
        if (status) then
            read(string, *, err=9999) moon_diampts
            if (moon_diampts .lt. 0.) goto 9999
        else
            moon_diampts = 0.d0
        end if

c***********************************************************************
c Determine whether disks should be blank
c***********************************************************************

        keyword = 'BLANK_DISKS'
        status = WWW_GetKey('blank', string)
        if (string .eq. 'Yes') then
            blank_disks = .TRUE.
        else if (string .eq. 'No' .or. string .eq. '') then
            blank_disks = .FALSE.
        else
            goto 9999
        end if

c***********************************************************************
c Determine whether meridians should be highlighted
c***********************************************************************

        keyword = 'HIGHLIGHTED_MERIDIANS'
        status = WWW_GetKey('meridians', string)
        if (string .eq. 'Yes') then
            meridian_pts = 1.3d0
        else if (string .eq. 'No' .or. string .eq. '') then
            meridian_pts = 0.d0
        else
            goto 9999
        end if

c***********************************************************************
c Generate output table of body positions
c***********************************************************************

c Locate planet
        call RSPK_BodyRaDec(obs_time, PLANET_ID, planet_ra, planet_dec)
        cosdec = cos(planet_dec)

c Write RA/dec table
        write(*,10) ' '
        write(*,10) 'Field of View Description (J2000)'
        write(*,10) '---------------------------------'
        write(*,10) ' '
        if (obs_isc .eq. 0) then
            write(*,11)
11          format(1x, 4x, 'Body', 10x, 'RA', 18x, 'Dec',
     &             18x, 'RA (deg)', 3x, 'Dec (deg)',
     &             4x, 'dRA (")', 3x, 'dDec (")')
        else
            write(*,12)
12          format(1x, 4x, 'Body', 10x, 'RA', 18x, 'Dec',
     &          18x, 'RA (deg)', 3x, 'Dec (deg)',
     &          3x, 'dRA (deg)', 2x, 'dDec (deg)')
        end if

        do 400 i = 0, NMOONS
            if (.not. moon_flags(i)) goto 400

            call RSPK_BodyRaDec(obs_time, MOON_IDS(i), ra, dec)
            if (obs_isc .eq. 0) then
                ddec = DPR() * 3600.d0 * (dec - planet_dec)
                dra  = DPR() * 3600.d0 * (ra  - planet_ra)
                if (dra .lt. -0.5*MAXSECS) dra = dra + MAXSECS
                if (dra .gt.  0.5*MAXSECS) dra = dra - MAXSECS
            else
                ddec = DPR() * (dec - planet_dec)
                dra  = DPR() * (ra  - planet_ra)
                if (dra .lt. -180.) dra = dra + 360.
                if (dra .gt.  180.) dra = dra - 360.
            end if
            dra = dra * cosdec

            ra = ra * DPR()
            dec = dec * DPR()
            call DMS_string(ra/15.d0,  'hms', 4, ra_string)
            call DMS_string(dec, 'dms', 3, dec_string)

            if (obs_isc .eq. 0) then
                write(*,13,err=400) MOON_IDS(i), MOON_NAMES(i),
     &                              ra_string, dec_string,
     &                              ra, dec, dra, ddec
13              format(1x, i3, 1x, a10, 3x, a18, 2x, a18,
     &                 2x, f10.6, f12.6, 1x, f10.4, 1x, f10.4)
            else
                write(*,14,err=400) MOON_IDS(i), MOON_NAMES(i),
     &                              ra_string, dec_string,
     &                              ra, dec, dra, ddec
14              format(1x, i3, 1x, a10, 3x, a18, 2x, a18,
     &                 2x, f10.6, f12.6, 1x, f11.6, 1x, f11.6)
            end if
400     continue

c Write body geometry table
        write(*,10) ' '
        write(*,22)
        write(*,23) LONDIR, LONDIR
22      format(19x,'Sub-Observer  ',7x,'Sub-Solar     ')
23      format(5x,'Body',10x,
     &         2('Lon(deg',a1,') Lat(deg)',3x),
     &         'Phase(deg)',3x,'Distance(10^6 km)')
        do 401 i = 0, NMOONS
            if (.not. moon_flags(i)) goto 401

            call RSPK_BodyLonLat(obs_time, MOON_IDS(i),
     &                           subobs_lon, subobs_lat,
     &                           subsol_lon, subsol_lat)
            call RSPK_BodyRanges(obs_time, MOON_IDS(i), sun_dist,
     &                                                  obs_dist)
            call RSPK_BodyPhase(obs_time, MOON_IDS(i), phase)

            write(*,24) MOON_IDS(i), MOON_NAMES(i),
     &                      subobs_lon * DPR(), subobs_lat * DPR(),
     &                      subsol_lon * DPR(), subsol_lat * DPR(),
     &                      phase * DPR(), obs_dist/1.d6
401     continue
24      format(1x, i3, 1x, a10,
     &         1x,f10.3,f10.3,1x,f10.3,f10.3,f13.5,f15.6)

c Write ring geometry values
        write(*,10) ' '

        call RSPK_RingOpen(obs_time, obs_b, sun_b, sun_db, isdark,
     &          obs_long, sun_long)
        call RSPK_Phase(obs_time, phase)
        call RSPK_Ranges(obs_time, sun_dist, obs_dist)
        litstr = '(lit)'
        if (isdark) litstr = '(unlit)'

        write(*,31) '  Ring sub-solar latitude (deg): ',
     &          DPR() * sun_b,
     &          DPR() * (sun_b-sun_db),
     &          DPR() * (sun_b+sun_db)
        write(*,32) ' Ring plane opening angle (deg): ',
     &          DPR() * obs_b,
     &          litstr(1:LASTNB(litstr))
        write(*,32) '  Ring center phase angle (deg): ',
     &          phase * DPR()
        write(*,32) '      Sub-solar longitude (deg): ',
     &          DPR() * sun_long,
     &          'from ring plane ascending node'
        write(*,32) '   Sub-observer longitude (deg): ',
     &          DPR() * obs_long
        write(*,*)
        write(*,32) '       Sun-planet distance (AU): ',
     &          sun_dist/AU
        write(*,32) '  Observer-planet distance (AU): ',
     &          obs_dist/AU
        write(*,33) '       Sun-planet distance (km): ',
     &          sun_dist/1.d6,
     &          'x 10^6'
        write(*,33) '  Observer-planet distance (km): ',
     &          obs_dist/1.d6,
     &          'x 10^6'
        write(*,33) '        Light travel time (sec): ',
     &          obs_dist/CLIGHT(0.d0)
        write(*,*)
31      format(1x, a, f9.5, '  (', f9.5, '  to ', f9.5, ')')
32      format(1x, a, f9.5, 2x, a)
33      format(1x, a, f12.6, 1x, a)

c***********************************************************************
c Determine ring offsets
c***********************************************************************

c Locate barycenter
        call RSPK_ObsLoc(obs_time, obs_pv)
        call SPKAPP(9, obs_time, 'J2000', obs_pv, 'LT',
     &              bary_dpv, dt)
        call SPKAPP(999, obs_time, 'J2000', obs_pv, 'LT',
     &              planet_dpv, dt)

        ring_offsets(1,1) = 0.
        ring_offsets(2,1) = 0.
        ring_offsets(3,1) = 0.

        do i = 2, NRINGS
            ring_offsets(1,i) = bary_dpv(1) - planet_dpv(1)
            ring_offsets(2,i) = bary_dpv(2) - planet_dpv(2)
            ring_offsets(3,i) = bary_dpv(3) - planet_dpv(3)
        end do

c***********************************************************************
c Set up plot title and captions
c***********************************************************************

        status = WWW_GetKeyUnsanitized('title', title)

        ncaptions = 1
        lcaptions(ncaptions) = 'Time (UTC):'
        status = WWW_GetKey('time', rcaptions(ncaptions))

        ncaptions = ncaptions + 1
        lcaptions(ncaptions) = 'Ephemeris:'
        status = WWW_GetKey('ephem', rcaptions(ncaptions))
        rcaptions(ncaptions) = rcaptions(ncaptions)(5:)

        ncaptions = ncaptions + 1
        lcaptions(ncaptions) = 'Viewpoint:'
        status = WWW_GetKey('viewpoint', string)
        if (string .eq. 'latlon') then
            string = '('
            status = WWW_GetKey('latitude', string(2:))
            string(LASTNB(string)+1:) = ','
            status = WWW_GetKey('longitude', string(LASTNB(string)+2:))
            status = WWW_GetKey('lon_dir',   string(LASTNB(string)+2:))
            string(LASTNB(string)+1:) = ','
            status = WWW_GetKey('altitude', string(LASTNB(string)+2:))
            string(LASTNB(string)+1:) = ')'
            rcaptions(3) = string
        else
            status = WWW_GetKey('observatory', string)
            status = WWW_GetKey('sc_trajectory', string2)
            if (status) then
                string(LASTNB(string)+1:) =
     &                  ' (' // string2(5:LASTNB(string2)) // ')'
            end if
            rcaptions(ncaptions) = string
        end if

        ncaptions = ncaptions + 1
        lcaptions(ncaptions) = 'Moon selection:'
        status = WWW_GetKey('moons', string)
        rcaptions(ncaptions) = string(5:)

        status = WWW_GetKey('moremoons', string2)
        if (status) then
            rcaptions(ncaptions) = string(5:LASTNB(string)) //
     &                             ' + ' // string2
        end if

        ncaptions = ncaptions + 1
        lcaptions(ncaptions) = 'Ring selection:'
        status = WWW_GetKey('rings', rcaptions(ncaptions))

        ncaptions = ncaptions + 1
        call RSPK_BodyLonLat(obs_time, center_body_id,
     &                       subobs_lon, subobs_lat,
     &                       subsol_lon, subsol_lat)
        lcaptions(ncaptions) =
     &          center_body_name(1:LASTNB(center_body_name)) //
     &          ' center (lon,lat):'
        write(rcaptions(ncaptions),80) subobs_lon * DPR(), LONDIR,
     &                                 subobs_lat * DPR()
80      format('(',f7.3,'\260 ',a1,',',f7.3,'\260)')

        ncaptions = ncaptions + 1
        call RSPK_BodyPhase(obs_time, center_body_id, phase)
        lcaptions(ncaptions) =
     &          center_body_name(1:LASTNB(center_body_name)) //
     &          ' phase angle: '
        write(rcaptions(ncaptions),81) phase * DPR()
81      format(f7.3,'\260')

c***********************************************************************
c Write Postscript file
c***********************************************************************

c Get Postscript filename
        call WWW_GetEnv('VIEWER_POSTFILE', filename)

c Planet is not a moon for this purpose
        moon_flags(0) = .FALSE.

        call RSPK_DrawView(obs_time, fov, center_ra, center_dec,
     &                  blank_disks, meridian_pts,
     &          NMOONS+1, moon_flags, MOON_IDS, MOON_NAMES,
     &                  moon_labelpts, moon_diampts,
     &          NRINGS, ring_flags, ring_rads, ring_elevs,
     &                  ring_eccs, ring_incs, ring_peris, ring_nodes,
     &                  ring_offsets,
     &                  ring_opaqs, ring_dashed, ring_method,
     &          NARCS, arc_flags, arc_rings, arc_minlons, arc_minlons,
     &                  arc_widthpts,
     &          nstars, star_ras, star_decs, star_names, star_labels,
     &                  STAR_DIAMPTS,
     &          title, ncaptions, lcaptions, rcaptions, 108.d0,
     &          filename)

        call exit

c***********************************************************************
c Handle invalid input
c***********************************************************************

9999    continue
        write(*,90) 'Invalid value found for variable ',
     &                  keyword(1:LASTNB(keyword)),
     &                  string(1:LASTNB(string))
90      format(1x, a, a, ': ', a)
        call exit

        end

c***********************************************************************
