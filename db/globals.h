// Global types and constants for EDASkel, a sample EDA app
// Copyright (C) 2010 Jeffrey Elliot Trull <edaskel@att.net>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

// Some definitions I'll need "everywhere" if only to avoid passing around strings
// This poses an interesting conundrum, though: how do I keep an enum generic enough
// to keep the templated database/library classes compatible with the parser, etc.?

#ifndef EDASKEL_GLOBALS_H
#define EDASKEL_GLOBALS_H

namespace EDASkel {

enum SiteSymmetry { SITESYM_X, SITESYM_Y, SITESYM_R90 };
enum SiteClass { SITECLASS_COVER, SITECLASS_RING, SITECLASS_BLOCK,
		 SITECLASS_PAD, SITECLASS_PAD_INPUT, SITECLASS_PAD_OUTPUT, SITECLASS_PAD_INOUT,
		 SITECLASS_PAD_POWER, SITECLASS_PAD_SPACER,
		 SITECLASS_CORE, SITECLASS_CORE_FEEDTHRU, SITECLASS_CORE_TIEHIGH, SITECLASS_CORE_TIELOW,
		 SITECLASS_ENDCAP_PRE, SITECLASS_ENDCAP_POST, SITECLASS_ENDCAP_TOPLEFT,
		 SITECLASS_ENDCAP_TOPRIGHT, SITECLASS_ENDCAP_BOTTOMLEFT, SITECLASS_ENDCAP_BOTTOMRIGHT };

 enum Orient { ORIENT_N, ORIENT_E, ORIENT_W, ORIENT_S,
	       ORIENT_FN, ORIENT_FE, ORIENT_FW, ORIENT_FS };
}
#endif
